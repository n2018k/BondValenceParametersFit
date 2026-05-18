import argparse
import os
from pathlib import Path

from bond_valence_processor import BondValenceProcessor


def _normalize_cations(*, cations: list[str], cation: list[str]) -> list[str]:
    tokens: list[str] = []
    tokens.extend(cations or [])
    tokens.extend(cation or [])

    normalized: list[str] = []
    seen: set[str] = set()
    for token in tokens:
        for part in str(token).split(","):
            sym = part.strip()
            if not sym:
                continue
            if sym in seen:
                continue
            seen.add(sym)
            normalized.append(sym)
    return normalized


def main() -> None:
    # Ensure output goes to this repo's ./res even if invoked from another directory.
    os.chdir(Path(__file__).resolve().parent)

    parser = argparse.ArgumentParser(description="Bond valence parameter fitting (batch or single MP material).")
    parser.add_argument("--api-key", default=os.environ.get("MP_API_KEY", ""), help="Materials Project API key (or set MP_API_KEY).")
    parser.add_argument(
        "--cations",
        nargs="+",
        default=[],
        help="One or more cation element symbols (e.g., Li Na K or Li,Na,K).",
    )
    parser.add_argument(
        "--cation",
        action="append",
        default=[],
        help="Deprecated alias for --cations (can be provided multiple times).",
    )
    parser.add_argument("--anion", default="O", help="Anion element symbol (e.g., O).")
    parser.add_argument(
        "--mp-id",
        default=None,
        help="Optional Materials Project material ID (e.g., mp-123). If set, only this structure is processed.",
    )
    parser.add_argument(
        "--algos",
        nargs="+",
        default=["shgo", "brute", "diff", "dual_annealing", "direct"],
        help="Optimization algorithms to run.",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Write debug grids to Excel instead of fitting R0/B.",
    )
    parser.add_argument(
        "--enforce-constraints",
        action="store_true",
        help="Enforce physical constraints on fitted parameters (no Sij tweaking unless --allow-sij-tweaks is set).",
    )
    parser.add_argument(
        "--allow-sij-tweaks",
        action="store_true",
        help="If constraint checks fail, retry by tweaking Sij on the shortest target bond (implies --enforce-constraints).",
    )
    args = parser.parse_args()

    if not args.api_key:
        raise SystemExit("Missing API key. Pass --api-key or set MP_API_KEY.")

    cations = _normalize_cations(cations=args.cations, cation=args.cation)
    if not cations:
        raise SystemExit("No cations provided. Pass --cations (or the deprecated --cation).")

    enforce_constraints = bool(args.enforce_constraints or args.allow_sij_tweaks)
    processor = BondValenceProcessor(
        args.api_key,
        args.algos,
        cations,
        args.anion,
        enforce_constraints=enforce_constraints,
        allow_sij_tweaks=bool(args.allow_sij_tweaks),
    )
    for cation in cations:
        processor.process_cation_system(cation, args.anion, mp_id=args.mp_id, debug=args.debug)


if __name__ == "__main__":
    main()
