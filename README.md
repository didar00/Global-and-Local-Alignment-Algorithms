### Global and Local Alignment Algorithms for Biomolecular Sequences
Implementations of Needleman-Wunsch (global alignment) and Smith-Waterman (local alignment) algorithms using Python.

Command line:
```
# G for global aligment
# L for local alignment
python main.py --input=input.txt --gap_extension=-5 --gap_opening=-10 --score_matrix=scoring_matrices/BLOSUM62.txt --alg=G
```
