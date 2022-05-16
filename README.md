# ChemosGeneAnnotation

Set of scripts to identify chemsensory related genes in unannotated genomes. Using raw data from 55 lepidopteran and dipteran species, this pipeline can be run iteratively for annotation and confirmation of genes involved in chemosensation in insects. It uses previously identified protein sequences for odorant receptors (ORs), gustatory receptors (GRs), ionotropic receptors (IRs), odorant-binding proteins (OBPs), chemosensory proteins (CSPs) and sensory neuron membrane proteins (SNMPs).

## Requirements

* python3
* BLAST+ suite - [install for command line](https://www.ncbi.nlm.nih.gov/books/NBK279671/)
* Exonerate - [install](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)
* Pfam for command line

## Usage

ChemosGeneAnnotation is a collection of python scripts that is currently run with a string of commands.

To run, download this repo locally using `git clone https://github.com/PeterMulhair/HbxFinder.git`

### Step 1: 