#!/usr/bin/python3
from subprocess import call as unix
from Bio import SeqIO
from joblib import Parallel, delayed
from collections import defaultdict
import argparse
import inflect
import os
import glob

'''
Run Exonerate on extracted contig regions
from BLAST analysis. Annotates full gene
with exon-intron boundaries for each of 
the CS genes.
'''

parse = argparse.ArgumentParser()

parse.add_argument("--input",type=str, help="path to input chemosensory genes to use in search",required=True)

args = parse.parse_args()


gene_list = ['CSP', 'GOBP', 'GR', 'IR', 'SNMP', 'OBP', 'OR']       
    
def exonerate(taxa, rounds):

    if rounds == 'NULL':
        output_path = 'results/exonerate_output'
    else:
        output_path = 'results/exonerate_output_' + rounds
        
    CS_gene = taxa.split('/')[-1].split('.')[0].split('_')[-1]
    if CS_gene in gene_list:

        gene_fasta = glob.glob(args.input + '*_' + CS_gene + '*')

        sp_out = taxa.split('/')[-1].split('.')[0].split('_')[0]
        gene_sp_out = sp_out + '_' + CS_gene + '.gff'
        
        if len(gene_fasta) > 0:
            unix('exonerate --model protein2genome --forcegtag --percent 10 --maxintron 8000 ' + gene_fasta[0] + ' ' + taxa + ' --showtargetgff yes --bestn 1 --verbose 0 --showalignment no --showvulgar no >> ' + output_path + '/' + gene_sp_out, shell=True)
    


if os.path.isdir('results/blastout/round_two'):
    previous_runs = len(next(os.walk('results/blastout/'))[1])
    new_round = previous_runs
    p = inflect.engine()
    new_round_word = p.number_to_words(new_round)
    print('Running exonerate round', new_round_word)
    os.mkdir('results/exonerate_output_' + new_round_word)
    os.mkdir('results/exonerate_output_' + new_round_word + '/sequences')
    os.mkdir('results/exonerate_output_' + new_round_word + '/sequences/separated_genes')
    unix('cp src/exonerate_parser.py results/exonerate_output_' + new_round_word + '/', shell=True)
    unix('cp src/edit_gff.py results/exonerate_output_' + new_round_word + '/sequences', shell=True)
    unix('sed -i "s:/blastout/sp_blast_extract:/blastout/round_' + new_round_word + '/sp_blast_extract:g" results/exonerate_output_' + new_round_word + '/sequences/edit_gff.py', shell=True)
    unix('cp src/AvrilGffUtils.pm results/exonerate_output_' + new_round_word + '/sequences', shell=True)
    unix('cp src/convert_exonerate_gff_to_std_gff.pl results/exonerate_output_' + new_round_word + '/sequences', shell=True)
    unix('cp src/sep_translate_genes.py results/exonerate_output_' + new_round_word + '/sequences/separated_genes', shell=True)
    unix('cp src/gene_counts.py results/exonerate_output_' + new_round_word + '/sequences/separated_genes', shell=True)
    sp_chemo_blast = []
    for sp in glob.glob("results/blastout/round_" + new_round_word + "/sp_blast_extract/*.fasta"):
        sp_chemo_blast.append(sp)

    Parallel(n_jobs=10)(delayed(exonerate)(taxa, new_round_word) for taxa in sp_chemo_blast)
        
else:
    new_round_word = 'NULL'
    print('Running first exonerate run')
    os.mkdir('results/exonerate_output')
    sp_chemo_blast = []
    for sp in glob.glob("results/blastout/sp_blast_extract/*.fasta"):
        sp_chemo_blast.append(sp)

    Parallel(n_jobs=10)(delayed(exonerate)(taxa, new_round_word) for taxa in sp_chemo_blast)
    
