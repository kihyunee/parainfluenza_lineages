import argparse
import os



def parse_seqs_from_fasta(fasta_file):
    list_seqid = []
    list_sequence = []
    fr = open(fasta_file, 'r')
    line = fr.readline()
    while line != '':
        if line.strip().startswith('>'):
            seqid = line.strip()[1:].split(' ')[0]
            seq_l = []
            line = fr.readline()
            while line != '':
                seq_l.append(line.strip())
                line = fr.readline()
                if line.strip().startswith('>'):
                    break
            list_seqid.append(seqid)
            list_sequence.append(''.join(seq_l))
        else:
            line = fr.readline()
    fr.close()
    return [list_seqid, list_sequence]


dict_valid_ntbase = {}
dict_valid_ntbase['A'] = True;  dict_valid_ntbase['a'] = True
dict_valid_ntbase['T'] = True;  dict_valid_ntbase['t'] = True
dict_valid_ntbase['C'] = True;  dict_valid_ntbase['c'] = True
dict_valid_ntbase['G'] = True;  dict_valid_ntbase['g'] = True
dict_valid_ntbase['U'] = True;  dict_valid_ntbase['u'] = True


# input structure
parser = argparse.ArgumentParser(description = "Input specifications ")
parser.add_argument("--refseq-fasta", dest="refseq_fasta_file", required=True, type=str, help="genomic fasta file of the reference genome used to align")
parser.add_argument("--query-aln-fofn", dest="query_pwaln_fofn", required=True, type=str, help="two column fofn: col 1 = query genome accession; col 2 = path to the ref-query pairwise alignment fasta file")
parser.add_argument("--cov", dest="site_coverage_cutoff", required=True, type=float, help="[0 .. 1] alignment coverage cutoff to include the site")
parser.add_argument("--out", dest="output_msa_fasta", required=True, type=str, help="output multiple sequence alignment fasta")



args = parser.parse_args()
refseq_fasta_file = args.refseq_fasta_file
query_pwaln_fofn = args.query_pwaln_fofn
site_coverage_cutoff = args.site_coverage_cutoff
output_msa_fasta = args.output_msa_fasta


# query fofn
list_query_seqid = []
list_query_pwaln_file = []
fr = open(query_pwaln_fofn, 'r')
for line in fr:
    fields = line.strip().split("\t")
    list_query_seqid.append(fields[0])
    list_query_pwaln_file.append(fields[1])
fr.close()
num_query = len(list_query_seqid)
print("number of query sequences = " + str(num_query))


# reference sequence length
ref_fasta_seqid = []
ref_fasta_sequence = []
[ref_fasta_seqid, ref_fasta_sequence] = parse_seqs_from_fasta(refseq_fasta_file)
ref_seq_length = len(ref_fasta_sequence[0])
print("reference sequence length = " + str(ref_seq_length))

# coverage count slots
ref_site_cov_slot = [0]*ref_seq_length


# count up the coverage per ref site
print("count coverage per reference site")
for qidx in range(num_query):
    query_seqid = list_query_seqid[qidx]
    pwaln_file = list_query_pwaln_file[qidx]
    [aligned_seqids, aligned_sequences] = parse_seqs_from_fasta(pwaln_file)

    aligned_query_sequence = ''
    aligned_ref_sequence = ''
    for psidx in range(2):
        if aligned_seqids[psidx] == query_seqid:
            aligned_query_sequence = aligned_sequences[psidx]
        else:
            aligned_ref_sequence = aligned_sequences[psidx]
    
    num_covered_in_this_query = 0
    alignment_length = len(aligned_ref_sequence)
    ref_site_index = -1
    for aln_site_idx in range(alignment_length):
        if aligned_ref_sequence[aln_site_idx] != '-':
            ref_site_index += 1
            if aligned_query_sequence[aln_site_idx] in dict_valid_ntbase:
                ref_site_cov_slot[ref_site_index] += 1
                num_covered_in_this_query += 1

    print("  ...  " + query_seqid + " : alignment length = " + str(alignment_length) + " / covered reference site = " + str(num_covered_in_this_query))


# determine the valid site to include in output msa
print("determine the conserved sites")
ref_site_pass = [False]*ref_seq_length
num_site_pass = 0
for ref_site_index in range(ref_seq_length):
    count_cov = ref_site_cov_slot[ref_site_index]
    cov_rate = float(count_cov)/float(num_query)
    if cov_rate >= site_coverage_cutoff:
        ref_site_pass[ref_site_index] = True
        num_site_pass += 1
print("number of reference sites passing the coverage cutoff = " + str(num_site_pass))


# write up the msa
fw = open(output_msa_fasta, 'w')
for qidx in range(num_query):
    query_seqid = list_query_seqid[qidx]
    pwaln_file = list_query_pwaln_file[qidx]
    [aligned_seqids, aligned_sequences] = parse_seqs_from_fasta(pwaln_file)

    aligned_query_sequence = ''
    aligned_ref_sequence = ''
    for psidx in range(2):
        if aligned_seqids[psidx] == query_seqid:
            aligned_query_sequence = aligned_sequences[psidx]
        else:
            aligned_ref_sequence = aligned_sequences[psidx]
    fw.write(">" + query_seqid + "\n")

    alignment_length = len(aligned_ref_sequence)
    ref_site_index = -1
    for aln_site_idx in range(alignment_length):
        if aligned_ref_sequence[aln_site_idx] != '-':
            ref_site_index += 1
            if aligned_query_sequence[aln_site_idx] in dict_valid_ntbase:
                fw.write(aligned_query_sequence[aln_site_idx])
            else:
                fw.write("-")

    fw.write("\n")
fw.close()
