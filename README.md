# About the page

This is a supplementary material for a research article __Phylogenetic Lineage Dynamics of Global Parainfluenza Virus Type 3 Post-COVID-19 Pandemic__,
presented by Kihyun Lee, Kuenyoul Park, Heungsup Sung and Mi-Na Kim, unpublished yet. (We will update the article status information once the article go through peer review).

In this repository we provide the sequence data and the accompanying metadata files, plus the commands or scripts, that can be used to replicate the phylogenetic analyses provided in the article.

## Data files
*All data files can be found in the files branch of this repository or accessed through the link provided below, next to each item's title.*

(1) The tsv file that contains accession number and metadata (parsed from the genbank records) of all PIV3 sequences, regardless of the sequenced region (single gene - whole genome). The entries in this table include the sequences that were not analyzed in the subsequent steps.
* __PIV3_NCBI_accessions_metadata_genome_coverage.tsv__ --> [File]()

(2) Sequence fasta and the accession information for the subset of sequences that we considered as "whole genomes", i.e., covering all six known CDS regions of PIV-3 genome. For the time-scaled phylogenetic analysis, we used the items whose sampling date were resolved at least at month level.
* __list_analyzed_WG.20230423.tsv :__ [File]()
* __subset.month_resolved.WG.metadata_pseudoday.tsv :__ [File]()

(3) Sequence fasta and the accession information for the subset of sequences that we considered to cover the HN gene with >90% coverage. For the time-scaled phylogenetic analysis, we used the items whose sampling date were resolved at least at month level.
* __list_analyzed_HN.20230423.tsv :__ [File]()
* __subset.month_resolved.HN.metadata_pseudoday.tsv :__ [File]()
