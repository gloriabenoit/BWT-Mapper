# Short read aligner based on Burrows–Wheeler transform

October 2024 (M2 BI)

## Introduction
This project aims to recreate a simple short read mapper using Burrows–Wheeler transform. 
Using a single read and a fasta file, you can get all starting positions which correspond to your read in your sequence.

## Usage (command line interface)

```bash
python src/mapper.py sequence.fasta short_read
```

## Example

As an example, we have searched ``` ATGCAG ``` in [SARS-CoV-2 full genome](https://www.ncbi.nlm.nih.gov/nuccore/MN908947). We found 13 matches.   

```bash
python src/mapper.py data/WH1.fst ATGCAG
```  