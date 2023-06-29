# Synonymous and Nonsynonymous Mutation comparison

This simple workflow constructs RSV-A and RSV-B synonymous and nonsynonymous mutation graphs.
These graphs are scaled by gene length, tree length and multiplied based on the proportion of loci
 at which they may occur (3 for synonymous, as they can occur at every third codon, and 3/2 for nonsynonymous).

## Running the workflow

This workflow uses Snakemake. To run from the command line, run "snakemake --cores all" from this directory.
RSV-A and RSV-B can be specified as "a" or "b" in the Snakefile rule all input. 

## Workflow inputs

This workflow requires as input in the data folder:

* a nwk tree file to calculate total tree length

* a nucleotide mutation file (json format)

* an amino acid mutation file (json format)
