# bayestraits-wrapper
Perl wrapper scripts for running the BayesTraits program with gene presence/absence data across multiple genomes

## Prerequisites
**Perl libraries:** `BioPerl`, `Getopt::Long`, `Sort::Naturally` and `List::Util`.

**BayesTraits** (http://www.evolution.rdg.ac.uk/BayesTraits.html) should be discoverable in `$PATH`, and the files `DEP.command` and `INDEP.command` need to be present in the `cwd`.

The `*.command` files are required for running `BayesTraits`; see the manual for more information.

### Help

Run programs without commands or with `-h` flag to see the list of options.

## bayesTraitsWrapper_ML.pl

A Perl wrapper around the BayesTraits program, running the ML version of the analysis across genome-wide gene presence / absence data.

Takes as input:
* A fasta-formatted matrix of gene presence absence, where each column of the alignment represents an orthologous group, gene presence is indicated by a `1` and gene absence with `0`.
* A traits file, formatted as `[GENOMENAME]\t[1|0]`, where `1` = trait presence, `0` = trait absence.
* A tree file, nexus formatted (ensure taxa names correspond across all files).

Matrix file format:
```
>taxa1
010111010001111...
>taxa2
010101101111001...
>taxa3
101011100011010...
etc.
```

Traits file format:
```
taxa1 1
taxa2 0
taxa3 1
etc.
```

## bayesTraitsWrapper_ML_NULL.pl

Generates a **null distribution** for testing significance of results from `bayesTraitsWrapper_ML.pl` by permuting the trait values and|or the gene presence / absence data, with respect to the phylogeny of the input taxa. This essentially breaks any associations between gene presence and trait presence that may exist in the data.

## Utilities
### orthogroups_to_binary_fasta.pl

Convert Orthgroups.txt file, e.g., from OrthoFinder, to fasta format required above. Run `orthogroups_to_binary_fasta.pl -h` to see help.

### add_taxa_translate_to_nexus.pl

BayesTraits requires a nexus tree with a 'translate' block instead of taxa names in the newick string itself (see the BT manual, and see the example trees that ship with the BT program). The utility script `add_taxa_translate_to_nexus.pl` will insert a translate block to a tree file. Run `add_taxa_translate_to_nexus.pl -h` to see options.

### get_significant_seqs_from_LR.pl

Takes the table file result of `bayesTraitsWrapper_ML.pl` and a user-input LRStat threshold and will returns representative sequences in fasta format for each site that shows an LR > LRStat threshold. Currently picks simply the **longest sequence** from any significantly associated OG as the representative. Also prints a 1/0 matrix table of just the selected sites. Run without commands to see help.

## References

* Pagel (1994) Proc R Soc Lond B 255:37
* Barker & Pagel (2005) PLoS Comp Biol 1:e3
* Nowell RW et al., (2016) Mol Plant Pathol. 17:1409.
* http://www.evolution.rdg.ac.uk/BayesTraits.html
