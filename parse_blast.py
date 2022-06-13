#!/usr/bin/python3
import os
import glob
import argparse
import inflect
from Bio import SeqIO
from collections import defaultdict
from joblib import Parallel, delayed

'''
Following TBLASTN search for chemosensory genes, 
this script parses the blast outputs to extract 
regions annotated with putative CS genes.The 
extracted sequences are then put through Exonerate 
later for annotation.
'''


parse = argparse.ArgumentParser()

parse.add_argument("--taxa",type=str, help="name of species blast result to parse",required=True)
parse.add_argument("--path",type=str, help="path to genomes in fasta format",required=True)

args = parse.parse_args()

output_list = []
for output in glob.glob(args.taxa + '*fa'):
    output_list.append(output)

num_sp = len(output_list)
print('Parsing BLAST output for', num_sp, 'species')

#For each CS gene used in blast, pull out the BLAST coordinates from the output files
gene_list = ['CSP', 'GOBP', 'GR', 'IR', 'SNMP', 'OBP', 'OR']

#for CS_gene in gene_list:
def pull_CSseq(CS_gene, rounds):
    for taxa in output_list:
        taxa_file = taxa.split('/')[-1]
        sp_out = taxa.split('/')[-1].split('_')[2].split('.')[0]

        assemb = taxa.split('/')[-1].split('.blast')[0]
        assembly = glob.glob(args.path + assemb + '.f*')
        contig_seq = {}
        with open(assembly[0]) as f:
            for record in SeqIO.parse(f, 'fasta'):
                contigID = record.id    
                seq = str(record.seq)
                contig_seq[contigID] = seq

                
        if rounds == 'NULL':
            new_file = 'results/blastout/sp_blast_extract/' + sp_out + '_' + CS_gene + '.fasta'
        else:
            new_file = 'results/blastout/round_' + rounds + '/sp_blast_extract/' + sp_out + '_' + CS_gene + '.fasta'

        #with open('results/blastout/sp_blast_extract/' + sp_out + '_' + CS_gene + '.fasta', 'w') as outF:
        with open(new_file, 'w') as outF:
            with open(taxa) as f:
                gene_contigs = defaultdict(list)
                contig_seq_hit = defaultdict(list)

                for line in f:
                    lines = line.split('\t')
                    #geneID = lines[0].split('|')[1]
                    gene = lines[0]#.split('0')[0]
                    if CS_gene in gene:

                        sp_contig = lines[1]
                        qstart = lines[5]
                        qend = lines[6]
                        qhit = int(qend) - int(qstart)
                        qlen = lines[7]
                        q_cov = qhit/int(qlen)
                        q_cov = 100*q_cov

                        s_start = lines[8]
                        s_end = lines[9]
                        s_len = lines[10]

                        contig_seq_hit[sp_contig].append(int(s_start))
                        contig_seq_hit[sp_contig].append(int(s_end))
                        gene_contigs[gene].append(sp_contig)


                for contig_hit, nuc_hit in contig_seq_hit.items():
                    nuc_start = min(nuc_hit) - 100
                    nuc_end = max(nuc_hit) + 100 #Regions within a contig that had multiple CS gene hits are taken out in a single large sequence for annotation

                    #print(contig_hit, nuc_start, nuc_end)

                    genome_contig = contig_seq[contig_hit]

                    blast_hit = []
                    for i in range(nuc_start, nuc_end):
                        try:
                            nuc = genome_contig[i]
                            blast_hit.append(nuc)
                        except:
                            continue
                    blast_hit = ''.join(blast_hit)

                    outF.write('>' + contig_hit + '\n' + blast_hit + '\n')

                
if os.path.isdir('results/blastout/sp_blast_extract'):
    previous_runs = len(next(os.walk('results/blastout/'))[1])
    new_round = previous_runs
    p = inflect.engine()
    new_round_word = p.number_to_words(new_round)
    os.makedirs('results/blastout/round_' + new_round_word + '/sp_blast_extract', exist_ok=True)
    Parallel(n_jobs=7)(delayed(pull_CSseq)(gene, new_round_word) for gene in gene_list)
else:
    new_round_word = 'NULL'
    os.makedirs('results/blastout/sp_blast_extract', exist_ok=True)
    Parallel(n_jobs=7)(delayed(pull_CSseq)(gene, new_round_word) for gene in gene_list)
