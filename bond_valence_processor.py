from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from tqdm import tqdm
import os
import numpy as np
import json
import math
import re
from mp_api.client import MPRester
from BVparams_search import TheoreticalBondValenceSolver, BVParamSolver
import pandas as pd

@dataclass
class MaterialData:
    material_id: str
    possible_species: List[str]
    structure_graph: Optional[object] = None
    formula_pretty: Optional[str] = None

class BondValenceProcessor:
    def __init__(
        self,
        api_key: str,
        algos: List[str],
        cations: List[str],
        anion: str,
        *,
        enforce_constraints: bool = False,
        allow_sij_tweaks: bool = False,
    ):
        self.api_key = api_key
        self.algos = algos
        self.cations = cations
        self.anion = anion
        self.enforce_constraints = bool(enforce_constraints)
        self.allow_sij_tweaks = bool(allow_sij_tweaks)
        self._ensure_directories()

    def _compute_debug_grid(
        self,
        *,
        cation: str,
        anion: str,
        bond_type_list: List[str],
        network_valence: Dict[str, float],
        bond_lengths: Dict[str, float],
        R0_min: float = 1.0,
        R0_max: float = 2.5,
        R0_step: float = 0.05,
        B_min: float = 0.1,
        B_max: float = 0.8,
        B_step: float = 0.05,
    ) -> pd.DataFrame:
        target_bonds = [
            e for e in bond_type_list
            if re.split(r'\d+', e)[0] == cation and re.split(r'\d+', e)[1] == anion
        ]
        if not target_bonds:
            return pd.DataFrame(columns=["R0", "B", "SSE"])

        sij = np.array([network_valence[b] for b in target_bonds], dtype=float)
        if np.any(sij <= 0):
            return pd.DataFrame(columns=["R0", "B", "SSE"])

        rij = np.array([bond_lengths[b] for b in target_bonds], dtype=float)
        log_sij = np.log(sij)

        n_R0 = int(round((R0_max - R0_min) / R0_step)) + 1
        n_B = int(round((B_max - B_min) / B_step)) + 1
        R0_values = R0_min + np.arange(n_R0, dtype=float) * R0_step
        B_values = B_min + np.arange(n_B, dtype=float) * B_step

        rows = []
        for B in B_values:
            base = (-B * log_sij) - rij
            for R0 in R0_values:
                e = R0 + base
                sse = float(np.sum(e * e))
                rows.append((float(R0), float(B), sse))

        return pd.DataFrame(rows, columns=["R0", "B", "SSE"])

    def _ensure_directories(self) -> None:
        """Ensure required directories exist for all cations"""
        Path("res").mkdir(exist_ok=True)
        for cation in self.cations:
            cation_dir = Path(f"res/{cation}{self.anion}")
            cation_dir.mkdir(exist_ok=True)
            (cation_dir / "params").mkdir(exist_ok=True)
            (cation_dir / "R0Bs").mkdir(exist_ok=True)
            (cation_dir / "no_solu").mkdir(exist_ok=True)
            for algo in self.algos:
                (cation_dir / "R0Bs" / algo).mkdir(exist_ok=True)

    def get_possible_species(
        self,
        save_dir: str,
        docs: List[MaterialData],
        *,
        include_missing: bool = False,
    ) -> List[str]:
        """Extract possible species from materials documents.

        If `include_missing` is True, material IDs are included even when
        `possible_species` is missing/empty (saved as an empty list), so downstream
        steps can still run using the fallback charge map.
        """
        species_data = {}
        for doc in tqdm(docs, desc='Getting possible species'):
            if doc.possible_species:
                species_data[doc.material_id] = doc.possible_species
            elif include_missing:
                species_data[doc.material_id] = []
        
        output_file = Path(save_dir) / "params" / "dict_matID_possible_species.json"
        with output_file.open('w') as f:
            json.dump(species_data, f)
            
        return list(species_data.keys())

    def process_cation_system(
        self,
        cation: str,
        anion: str,
        *,
        mp_id: Optional[str] = None,
        debug: bool = False,
    ) -> None:
        """Process either a full cation-anion system (batch) or a single MP structure.

        If `mp_id` is provided, only that Materials Project material is processed.
        """
        if mp_id:
            print(f'Processing {cation}-{anion} for {mp_id}...')
        else:
            print(f'Processing {cation}-{anion} system...')
        
        # Setup processing pipeline
        docs = self._download_materials_data(cation, anion, material_ids=[mp_id] if mp_id else None)
        res_dir = f'res/{cation}{anion}'
        mids = self.get_possible_species(res_dir, docs, include_missing=bool(mp_id))
        if not mids:
            print(f"No materials with possible species found for {cation}-{anion} system, skipping...")
            return
        bonds_docs = self._download_bonding_data([mp_id] if mp_id else mids)
        
        # Initialize data structures
        results = {
            'sij': {},
            'charges': {},
            'solved': set(),
            'no_solution_global': [],
            'no_solution_by_algo': {algo: [] for algo in self.algos},
        }
        
        # Load previous results if they exist
        results.update(self._load_previous_results(res_dir))
        
        # Process materials
        debug_excel_path = Path(res_dir) / "debug_grid.xlsx"
        if debug:
            with pd.ExcelWriter(debug_excel_path, engine="openpyxl") as writer:
                self._process_materials(
                    bonds_docs=bonds_docs,
                    results=results,
                    res_dir=res_dir,
                    cation=cation,
                    anion=anion,
                    debug=debug,
                    debug_excel_writer=writer,
                )
        else:
            self._process_materials(
                bonds_docs=bonds_docs,
                results=results,
                res_dir=res_dir,
                cation=cation,
                anion=anion,
                debug=debug,
                debug_excel_writer=None,
            )
        
        # Save final results
        self._save_results(res_dir, results['sij'], results['charges'])

    def _download_materials_data(
        self,
        cation: str,
        anion: str,
        material_ids: Optional[List[str]] = None,
    ) -> List[MaterialData]:
        """Download materials data from Materials Project.

        If `material_ids` is provided, fetch only those materials (ignoring the
        energy-above-hull filter).
        """
        with MPRester(api_key=self.api_key) as mpr:
            if material_ids:
                try:
                    return mpr.materials.summary.search(
                        material_ids=material_ids,
                        fields=['material_id', 'possible_species'],
                    )
                except TypeError:
                    # Some mp-api versions may not support `material_ids` in summary.search.
                    # Fall back to a broader query and filter client-side.
                    docs = mpr.materials.summary.search(
                        elements=[cation, anion],
                        fields=['material_id', 'possible_species'],
                    )
                    return [d for d in docs if getattr(d, "material_id", None) in set(material_ids)]
            return mpr.materials.summary.search(
                elements=[cation, anion],
                energy_above_hull=(0.000, 0.05),
                fields=['material_id', 'possible_species']
            )

    def _download_bonding_data(self, material_ids: List[str]) -> List[MaterialData]:
        """Download bonding data from Materials Project"""
        with MPRester(api_key=self.api_key) as mpr:
            return mpr.materials.bonds.search(
                material_ids=material_ids,
                fields=['material_id', 'structure_graph', 'formula_pretty']
            )

    def _load_previous_results(self, res_dir: str) -> Dict:
        """Load previously computed results"""
        results = {
            'solved': set(),
            'no_solution_global': [],
            'no_solution_by_algo': {algo: [] for algo in self.algos},
        }
        
        # Load solved materials and find common solved files across all algorithms
        solved_sets = []
        for alg in self.algos:
            alg_dir = Path(res_dir) / "R0Bs" / alg
            if alg_dir.exists():
                solved_files = {f.stem for f in alg_dir.glob("*.txt")}
                solved_sets.append(solved_files)
        
        # Find intersection of all solved files across algorithms
        if solved_sets:
            common_solved = set.intersection(*solved_sets)
            results['solved'].update(common_solved)
        
        # Load no-solution cases (per algorithm files).
        for alg in self.algos:
            no_solu_file = Path(res_dir) / "no_solu" / f"{alg}.txt"
            if not no_solu_file.exists():
                continue
            loaded = np.loadtxt(no_solu_file, dtype=str).tolist()
            if not loaded:
                continue
            if isinstance(loaded[0], str):
                loaded = [loaded]
            results['no_solution_by_algo'][alg].extend([tuple(row) for row in loaded])
            
        return results

    def _process_materials(
        self,
        bonds_docs: List[MaterialData],
        results: Dict,
        res_dir: str,
        cation: str,
        anion: str,
        *,
        debug: bool = False,
        debug_excel_writer: Optional[pd.ExcelWriter] = None,
    ) -> None:
        """Process each material in the dataset"""
        solver = TheoreticalBondValenceSolver(
            species_matID_path=str(Path(res_dir) / "params" / "dict_matID_possible_species.json")
        )
        
        for material in tqdm(bonds_docs, desc=f'Processing {cation}-{anion} materials'):
            if (not debug) and (material.material_id in results['solved']):
                continue
                
            # Compute Sij values
            sij_data = solver.get_sij(
                material.material_id,
                material.structure_graph.structure,
                material.structure_graph
            )
            
            # Store results
            results['sij'][material.material_id] = sij_data[0]
            results['charges'][material.material_id] = sij_data[3]

            if debug:
                if debug_excel_writer is None:
                    raise ValueError("debug_excel_writer is required when debug=True")
                network_valence, bond_types, bond_lengths, _ = sij_data
                if network_valence:
                    df = self._compute_debug_grid(
                        cation=cation,
                        anion=anion,
                        bond_type_list=bond_types,
                        network_valence=network_valence,
                        bond_lengths=bond_lengths,
                    )
                    df.to_excel(
                        debug_excel_writer,
                        sheet_name=str(material.material_id),
                        index=False,
                    )
                continue

            # Process with algorithms
            self._run_algorithms(
                sij_data=sij_data,
                material=material,
                results=results,
                res_dir=res_dir,
                cation=cation,
                anion=anion
            )

    def _run_algorithms(self, sij_data: Tuple, material: MaterialData, 
                       results: Dict, res_dir: str, cation: str, anion: str) -> None:
        """Run all algorithms on the material"""
        network_valence, bond_types, bond_lengths, _ = sij_data
        
        if not network_valence:
            no_solution_case = (
                material.material_id, cation, anion, 
                material.formula_pretty, 'no_network_sol'
            )
            results['no_solution_global'].append(no_solution_case)
            self._save_no_solution(res_dir, results['no_solution_global'], results['no_solution_by_algo'])
            return
            
        for algorithm in self.algos:
            solver = BVParamSolver(
                save_dir=res_dir,
                algo=algorithm,
                no_sol=results['no_solution_by_algo'][algorithm]
            )
            
            solved = solver.solve_R0Bs(
                cation=cation,
                anion=anion,
                bond_type_list=bond_types,
                networkValence_dict=network_valence,
                bondLen_dict=bond_lengths,
                materID=material.material_id,
                chem_formula=material.formula_pretty,
                R0_bounds=(0, 5),
                enforce_constraints=self.enforce_constraints,
                allow_sij_tweaks=self.allow_sij_tweaks,
                record_fit_log=True,
            )

            solution, fit_log = solved

            # If this is a single-equation target-bond case and we failed constraints,
            # record the full family of exact solutions using the ORIGINAL (untweaked) Sij.
            line_solution = None
            if not solution:
                target_bonds = [
                    e for e in bond_types
                    if re.split(r'\d+', e)[0] == cation and re.split(r'\d+', e)[1] == anion
                ]
                if len(target_bonds) == 1:
                    bond = target_bonds[0]
                    sij0 = network_valence.get(bond, None)
                    rij = bond_lengths.get(bond, None)
                    if sij0 is not None and rij is not None and sij0 > 0:
                        lnS = math.log(float(sij0))
                        line_solution = {
                            "bond": bond,
                            "Rij": float(rij),
                            "Sij": float(sij0),
                            "lnS": float(lnS),
                            "equation": "R0 = Rij + B*ln(Sij)",
                            "constraints": {
                                "B_min": 0.05,
                                "R0_ge_B_plus_eps": 0.1,
                            },
                        }

            # Persist attempt-by-attempt fit log (including residuals) for auditing.
            log_file = Path(res_dir) / "R0Bs" / algorithm / f"{material.material_id}.fitlog.json"
            with log_file.open("w") as f:
                json.dump(
                    {
                        "material_id": material.material_id,
                        "formula_pretty": material.formula_pretty,
                        "cation": cation,
                        "anion": anion,
                        "algorithm": algorithm,
                        "enforce_constraints": bool(self.enforce_constraints),
                        "allow_sij_tweaks": bool(self.allow_sij_tweaks),
                        "min_B": 0.05,
                        "eps_R0_ge_B": 0.1,
                        "sij_step_frac": 0.05,
                        "sij_max_frac": 0.20,
                        "attempts": fit_log,
                        "line_solution": line_solution,
                    },
                    f,
                    indent=2,
                )

            if solution:
                output_file = Path(res_dir) / "R0Bs" / algorithm / f"{material.material_id}.txt"
                np.savetxt(output_file, solution)
            else:
                self._save_no_solution(res_dir, results['no_solution_global'], results['no_solution_by_algo'])

    def _save_no_solution(self, res_dir: str, global_no_solution: List, no_solution_by_algo: Dict[str, List]) -> None:
        """Save no-solution cases per algorithm (plus global failures).

        Each `no_solu/<algo>.txt` contains:
        - all global failures (e.g., network solve failures)
        - failures encountered under that specific algorithm
        """
        global_rows = [tuple(row) for row in (global_no_solution or [])]
        global_set = set(global_rows)

        for algorithm in self.algos:
            algo_rows = [tuple(row) for row in (no_solution_by_algo.get(algorithm) or [])]
            merged = []
            seen = set()
            for row in global_rows + algo_rows:
                if row in seen:
                    continue
                seen.add(row)
                merged.append(row)

            output_file = Path(res_dir) / "no_solu" / f"{algorithm}.txt"
            if merged:
                np.savetxt(output_file, merged, fmt='%s')
            else:
                # Ensure the file exists but is empty.
                output_file.write_text("")

    def _save_results(self, res_dir: str, sij_data: Dict, charges: Dict) -> None:
        """Save final results to JSON files"""
        with open(Path(res_dir) / "dict_sijs.json", 'w') as f:
            json.dump(sij_data, f)
        
        with open(Path(res_dir) / "dict_charges.json", 'w') as f:
            json.dump(charges, f)


if __name__ == "__main__":
    # User-defined parameters
    user_cations = ['Li', 'Na', 'K', 'Rb', 'Cs']  # Can be modified by user
    user_anion = 'O'
    user_algos = ['shgo', 'brute', 'diff', 'dual_annealing', 'direct']
    api_key = "your_api_key"  # Should be provided by user
    
    processor = BondValenceProcessor(
        api_key=api_key,
        algos=user_algos,
        cations=user_cations,
        anion=user_anion,
    )
    
    for cation in user_cations:
        processor.process_cation_system(cation, user_anion)
