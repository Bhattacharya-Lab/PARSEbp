# PARSEbp: Pairwise Agreement-based RNA Scoring with Emphasis on Base Pairings

by Sumit Tarafder and Debswapna Bhattacharya

Codebase for our <ins>P</ins>airwise <ins>A</ins>greement-based <ins>R</ins>NA <ins>S</ins>coring with <ins>E</ins>mphasis on <ins>B</ins>ase <ins>P</ins>airings (PARSEbp).

## Installation
```
pip install PARSEbp
```

Or

```
git clone https://github.com/Bhattacharya-Lab/PARSEbp.git
cd PARSEbp
pip install .
```

Typical installation time should take less than a minute in a 64-bit Linux system.

## Usage

Instructions for running PARSEbp:

```python
# Import
from PARSEbp import parsebp

# Initialize
p = parsebp()

# Set target sequence to "" for sequence-agnostic scoring
p.set_target_sequnece("")

# Load a directory containing RNA 3D structures (.pdb files)
p.load_pdbs("Inputs/")

# Compute scores
score = p.score()

# Save the results
score.save("score.txt")
```

Additional functionality

```python
# Set scoring mode (default is 1)
p.set_mode(1)

# Set target sequence (only the models that exactly match the target sequence will be scored)
seq = "GGACACGAGUAACUCGUCUAUCUGCUGCAGGCUGCUUACGGUUUCGUCCGUGUUGCAGCCGAUCAUCAGAACAUCUAGGUUUCGUCCGGGUGUUACCGAAAGGUCAGAUGGAGAGCCUUGUCCC"
p.set_target_sequnece(seq)

# Set the number of threads for parallel computations (default is 50)
p.set_parallel_threads(50)

# Get score of a specific model
score.getScore("decoy_1.pdb")

# Get top-1 ranked model(s)
score.top1()

# Get top-N ranked decoys
score.topN(10)
```

Given a directory containing RNA 3D structures "Inputs" as input, PARSEbp predicts the quality score of each structure in the directory and saves the output in "Score.txt".

Score calculation for a typical RNA (~100 nucleotides) with ~200 3D structures takes ~30 seconds.

## Detailed Instructions

Follow the provided [notebook](./PARSEbp_colab.ipynb) for detailed explanation of the installation, scoring and score analysis.

## Datasets

- CASP16 3D decoy structures (all submitted predictions) are downloaded from [here](https://predictioncenter.org/download_area/CASP16/predictions/RNA/). 
- Ground truth scores for benchmarking PARSEbp is downloaded from [here](https://predictioncenter.org/casp16/results.cgi?tr_type=rna)

