# About the page

This is a supplementary material for a research article __Phylogenetic Lineage Dynamics of Global Parainfluenza Virus Type 3 Post-COVID-19 Pandemic__,
presented by Kihyun Lee, Kuenyoul Park, Heungsup Sung and Mi-Na Kim, unpublished yet. (We will update the article status information once the article go through peer review).

In this repository we provide the sequence data and the accompanying metadata files, plus the commands or scripts, that can be used to replicate the phylogenetic analyses provided in the article.

## Data files
*All data files can be found in the files branch of this repository or accessed through the link provided below, next to each item's title.*

(1) The tsv file that contains accession number and metadata (parsed from the genbank records) of all PIV3 sequences, regardless of the sequenced region (single gene - whole genome). The entries in this table include the sequences that were not analyzed in the subsequent steps.
* __PIV3_NCBI_accessions_metadata_genome_coverage.tsv__ --> [File](https://github.com/kihyunee/parainfluenza_lineages/blob/88605fda870c4ca9e3a51912fa13133941b3d091/PIV3_NCBI_accessions_metadata_genome_coverage.tsv)

(2) Sequence fasta and the accession information for the subset of sequences that we considered as "whole genomes", i.e., covering all six known CDS regions of PIV-3 genome. For the time-scaled phylogenetic analysis, we used the items whose sampling date were resolved at least at month level.
* __list_analyzed_WG.20230423.tsv__  --> [File](https://github.com/kihyunee/parainfluenza_lineages/blob/88605fda870c4ca9e3a51912fa13133941b3d091/list_analyzed_WG.20230423.tsv)
* __subset.month_resolved.WG.metadata_pseudoday.tsv__  --> [File](https://github.com/kihyunee/parainfluenza_lineages/blob/88605fda870c4ca9e3a51912fa13133941b3d091/subset.month_resolved.WG.metadata_pseudoday.tsv)
* __PIV3_nuc_20230423.HPIV3_ref_cds.blastn.all_gene_cov90.fasta__ --> [File](https://github.com/kihyunee/parainfluenza_lineages/blob/88605fda870c4ca9e3a51912fa13133941b3d091/PIV3_nuc_20230423.HPIV3_ref_cds.blastn.all_gene_cov90.fasta)

(3) Sequence fasta and the accession information for the subset of sequences that we considered to cover the HN gene with >90% coverage. For the time-scaled phylogenetic analysis, we used the items whose sampling date were resolved at least at month level.
* __list_analyzed_HN.20230423.tsv__  --> [File](https://github.com/kihyunee/parainfluenza_lineages/blob/88605fda870c4ca9e3a51912fa13133941b3d091/list_analyzed_HN.20230423.tsv)
* __subset.month_resolved.HN.metadata_pseudoday.tsv__  --> [File](https://github.com/kihyunee/parainfluenza_lineages/blob/88605fda870c4ca9e3a51912fa13133941b3d091/subset.month_resolved.HN.metadata_pseudoday.tsv)
* __PIV3_nuc_20230423.HN_gene.fasta__ --> [File](https://github.com/kihyunee/parainfluenza_lineages/blob/88605fda870c4ca9e3a51912fa13133941b3d091/PIV3_nuc_20230423.HN_gene.fasta)


## Scripts
(1) Whole genome multiple sequence alignment
We first aligned each whole genome sequence entry against a fixed reference genome. 
This was done by the following bash script, where we used *__ginsi__* command of MAFFT software to align each entry in the query fasta (containing multiple entries) vs. a fixed reference file (KF530231.fastsa).

```
#!/bin/bash

if [ ! -d ref_pw_ginsi ]; then
        mkdir ref_pw_ginsi
fi

query_wg_fasta=$1
n_query=$(grep -c ">" ${query_wg_fasta})
n_proc=0

while [ ${n_proc} -lt ${n_query} ]
do
        let n_proc=${n_proc}+1
        let qfrom=2*${n_proc}-1
        let qto=2*${n_proc}
        echo ""
        echo "###"
        echo "### ${n_proc}/${n_query} "
        echo "###    ${qfrom} .. ${qto}p"
        echo "###"
        echo ""
        sed -n ${qfrom},${qto}p ${query_wg_fasta} > tmp_query_seq.template
        query_acc="$(grep ">" tmp_query_seq.template | tr -d '>')"
        cat tmp_query_seq.template > ref_pw_ginsi/${query_acc}.template
        cat ../reference_ncbi_virus_genomes/HPIV3_ref_fasta/KF530231.fasta >> ref_pw_ginsi/${query_acc}.template
        ginsi ref_pw_ginsi/${query_acc}.template > ref_pw_ginsi/${query_acc}.aln
        rm tmp_query_seq.template

done
```

Then we prepared the file-of-file-names (fofn) tsv file indicating the path to ginsi output file (column 2) per input strain (column 1). The fofn file looks like this.
```
head ref_pw_ginsi.fofn
MK167032.2      ref_pw_ginsi/MK167032.2.aln
KY973582.2      ref_pw_ginsi/KY973582.2.aln
KY973581.2      ref_pw_ginsi/KY973581.2.aln
KY973579.2      ref_pw_ginsi/KY973579.2.aln
KY973577.2      ref_pw_ginsi/KY973577.2.aln
KY973575.2      ref_pw_ginsi/KY973575.2.aln
MF973173.1      ref_pw_ginsi/MF973173.1.aln
...
```

Then we joined the individual alignments into multiple alignment by a simple custom python script. The script can be found here, --> [python file](https://github.com/kihyunee/parainfluenza_lineages/blob/69a0340f2bbb9c12dacc25643534a2f269a0e6c0/pairwise_ref_aligns_to_msa.py)
```
python pairwise_ref_aligns_to_msa.py --refseq-fasta ../reference_ncbi_virus_genomes/HPIV3_ref_fasta/KF530231.fasta --query-aln-fofn ref_pw_ginsi.fofn --cov 0.9 --out ref_pw_ginsi.cov_90_sites_MSA.fasta
```

(2) Whole genome maximum likelihood phylogentic tree inference

(3) HN gene multiple sequence alignment

(4) HN gene maximum likelihood phylogenetic tree inference

(5) Codon-by-codon alignment per CDS
I used an in-house python script that takes nucleotide input sequences in fasta format and generate codon-aligned nucleotide alignment fasta file. 
This script is available in this repository's file section --> [python file](https://github.com/kihyunee/parainfluenza_lineages/blob/bb30d18b5bf9ce584a546855ed2d24090079029a/codon_alignment_from_cds_nt_fasta.py). 
You can also find the input CDS-by-CDS fasta files in the file section.
```
python codon_alignment_from_cds_nt_fasta.py --in genomic_gene_extract.452_genome.NP.fasta --out genomic_gene_extract.452_genome.NP.codon_align.aln --muscle [PATH_TO_MUSCLE5]
python codon_alignment_from_cds_nt_fasta.py --in genomic_gene_extract.452_genome.PCD.fasta --out genomic_gene_extract.452_genome.PCD.codon_align.aln --muscle [PATH_TO_MUSCLE5]
python codon_alignment_from_cds_nt_fasta.py --in genomic_gene_extract.452_genome.F.fasta --out genomic_gene_extract.452_genome.F.codon_align.aln --muscle [PATH_TO_MUSCLE5]
python codon_alignment_from_cds_nt_fasta.py --in genomic_gene_extract.452_genome.HN.fasta --out genomic_gene_extract.452_genome.HN.codon_align.aln --muscle [PATH_TO_MUSCLE5]
python codon_alignment_from_cds_nt_fasta.py --in genomic_gene_extract.452_genome.M.fasta --out genomic_gene_extract.452_genome.M.codon_align.aln --muscle [PATH_TO_MUSCLE5]
python codon_alignment_from_cds_nt_fasta.py --in genomic_gene_extract.452_genome.L.fasta --out genomic_gene_extract.452_genome.L.codon_align.aln --muscle [PATH_TO_MUSCLE5]
```

(6) SLAC on data monkey
You can easily replicate the site-by-site substitution rate analysis presented in the paper by going to [SLAC on Datamonkey](https://www.datamonkey.org/slac) (https://www.datamonkey.org/slac) and upload the codon-aligned fasta files created in the previous section.

(7) TempEst and BEAST on whole genome sequences

(8) TempEst and BEAST on HN gene sequences




