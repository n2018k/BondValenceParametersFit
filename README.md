# BondValenceParametersFit
A Python package for fitting bond valence parameters $R_0$ and $B$ for cation-anion pairs using Materials Project data. It includes two modules: (1) computing theoretical bond valence and (2) optimizing parameters by matching computed and empirical values. This tool refines bond valence analysis with a data-driven approach.

## How to fit BV parameters

```python
from bond_valence_processor import BondValenceProcessor

cations = ['Li'] # a list of cation species 
anion = 'O' # anion 
my_api_key = "your_api_key" # api key for materialsproject.org 
algos = ['shgo', 'brute', 'diff', 'dual_annealing', 'direct']
processor = BondValenceProcessor(my_api_key, algos, cations, anion)
    
for cation in cations:
    processor.process_cation_system(cation, anion)

## online bond valence parameter fitting web 
An online web interface for fitting bond valence parameters is available at Hugging Face: https://huggingface.co/spaces/nodameCL/BVSearchApp.
Users can upload a single CIF file to calculate the theoretical bond valence (S$_ij$) and obtain fitted bond valence parameters (R$_0$ and B) for target bonds with specified bond types. A snapshot of the interactive web application is displayed below. 
![interactive web app](figures/web_interface)    
```

