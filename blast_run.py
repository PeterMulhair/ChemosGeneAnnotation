#!/usr/bin/python3 
import os
import glob
import os.path
import inflect
import argparse
from joblib import Parallel, delayed
from subprocess import call as unix


parse = argparse.ArgumentParser()

parse.add_argument("--path",type=str, help="path to genomes in fasta format",required=True)
parse.add_argument("--threads",type=str, help="number of species to run in parallel eg the no. of genomes",required=True)
parse.add_argument("--input",type=str, help="raw input of chemosensory genes to use as search",required=True)

args = parse.parse_args()

os.makedirs('results', exist_ok=True)
os.makedirs('results/blastout', exist_ok=True)

db_list = []
sp_name_list = []
genome_list = glob.glob(args.path + "*.fasta") + glob.glob(args.path + "*.fa") + glob.glob(args.path + "*.fas") + glob.glob(args.path + "*.fna")
for fasta in genome_list:
        sp_fasta = fasta.split('/')[-1].split('.f')[0]
        sp_name_list.append(sp_fasta)
        db_list.append(fasta)

def run_blast(sp):
        print(sp)
        unix('makeblastdb -dbtype nucl -in ' + args.path + sp + '.f* -out ' + args.path + sp, shell=True)
        print('\nRunning BLAST search, this may take some time...\n')
        unix('tblastn -query raw/chemo_genes.fa -db ' + args.path + sp + ' -evalue 1e-5 -seg yes -max_target_seqs 5000 -num_threads 4 -outfmt "6 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen" -out results/blastout/' + sp + '.blastoutput.fa', shell=True)


def new_blast(sp, rounds):
        print(sp)
        print('\nRunning BLAST search, this may take some time...\n')
        unix('tblastn -query ' + args.input + ' -db ' + args.path + sp + ' -evalue 1e-5 -seg yes -max_target_seqs 5000 -num_threads 4 -outfmt "6 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen" -out results/blastout/round_' + rounds + '/' + sp + '.blastoutput.tsv', shell=True)



if os.path.isfile('results/blastout/' + sp_name_list[0] + '.blastoutput.fa'):
        previous_runs = len(next(os.walk('results/blastout/'))[1])
        new_round = previous_runs + 1
        p = inflect.engine()
        new_round_word = p.number_to_words(new_round)
        os.mkdir('results/blastout/round_' + new_round_word)
        print('New run: round', new_round_word)
        Parallel(n_jobs=int(args.threads))(delayed(new_blast)(sp_assem, new_round_word) for sp_assem in sp_name_list)
else:
        print('First run')
        Parallel(n_jobs=int(args.threads))(delayed(run_blast)(sp_assem) for sp_assem in sp_name_list)  
                
