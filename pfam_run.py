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
Use pfam database to annotate protein
domains in predicted genes. Code will
use last round of exonerate output as
input for domain annotation.
'''


gene_list = ['CSP', 'GOBP', 'GR', 'IR', 'SNMP', 'OBP', 'OR']       
gene_dict = {'CSP': ['OS-D'], 'GOBP': ['PBP_GOBP'], 'GR': ['Trehalose_recp', '7tm_7'], 'IR': ['Lig_chan', 'Lig_chan-Glu_bd'], 'SNMP': ['CD36'], 'OBP': ['PBP_GOBP'], 'OR': ['7tm_6']}

def pfam(taxa, rounds):
    output_path = 'results/exonerate_output_' + rounds + '/sequences/separated_genes/pfam_output/'
        
    CS_gene = taxa.split('/')[-1].split('.')[0].split('_')[-1]
    sp_gene = taxa.split('/')[-1].split('.')[0] + '_pfam.res'
    gene_domain = gene_dict[CS_gene]
    print(sp_gene)
    if CS_gene in gene_list:

        unix('perl ~/bin/PfamScan/PfamScan/pfam_scan.pl -fasta ' + taxa + ' -dir ~/bin/PfamScan/PfamScan/ -outfile ' + output_path + sp_gene, shell=True)
        unix('java -jar ~/bin/PfamScan/PfamScanner/PfamScanner.jar -in ' + output_path + sp_gene + ' -out ' + output_path + sp_gene + '.parsed -e 1e-3 -ne -p', shell=True)


        sp_gene_pass = []
        with open(output_path + sp_gene + '.parsed') as f:
            for record in SeqIO.parse(f, 'fasta'):
                gene = record.id
                domains = str(record.seq)
                domains = domains.split('\t')
                if any(gene_domains in domains for gene_domains in gene_domain):
                    sp_gene_pass.append(gene)

        with open(taxa) as f, open(output_path + taxa.split('/')[-1].split('.')[0] + '_pfam.fasta', 'w') as outF:
            for record in SeqIO.parse(f, 'fasta'):
                ID = record.id
                if ID in sp_gene_pass:
                    seq = str(record.seq)
                    outF.write('>' + ID + '\n' + seq + '\n')


if os.path.isdir('results/blastout/sp_blast_extract'):
    previous_runs = len(next(os.walk('results/blastout/'))[1])
    new_round = previous_runs
    p = inflect.engine()
    new_round_word = p.number_to_words(new_round)
    print('Running pfam on round', new_round_word, 'output\n')
    os.mkdir('results/exonerate_output_' + new_round_word + '/sequences/separated_genes/pfam_output')

    sp_chemo_exonerate = []
    for sp_gene in glob.glob('results/exonerate_output_' + new_round_word + '/sequences/separated_genes/*fasta'):
        sp_chemo_exonerate.append(sp_gene)

    Parallel(n_jobs=10)(delayed(pfam)(taxa, new_round_word) for taxa in sp_chemo_exonerate)

else:
    print('Cannot run pfam - need to run blast+exonerate first')
    sys.exit()
    
