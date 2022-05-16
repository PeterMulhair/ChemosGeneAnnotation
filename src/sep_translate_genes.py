#!/usr/bin/env python3
import glob
from joblib import Parallel, delayed
from subprocess import call as unix
from Bio import SeqIO, Align
import os

'''
Get non-overlapping annotated exonerate 
genes, separate into individual gene 
files and translate to amino acid sequences
retaining the largest peptide as representative
'''

aligner = Align.PairwiseAligner()

##Create dictionary of chemosensory gene IDs and AA sequences
chemo_gene_seq_dict = {}
for chemo_gene in glob.glob('/home/zoo/zool2500/repos/ChemosGeneAnnotation/raw/*fasta'):
    gene_name = chemo_gene.split('/')[-1].split('.')[0]
    seq_list = []
    if gene_name != 'chemo_genes':
        with open(chemo_gene) as f:
            for record in SeqIO.parse(f, 'fasta'):
                header = record.description
                seq = str(record.seq)
                seq_list.append(seq)
    longest_seq = max(seq_list, key=len)            
    chemo_gene_seq_dict[gene_name] = longest_seq

for fa in glob.glob("../*.fa"):
    sp = fa.split('/')[-1].split('_')[0]
    sp_gene = fa.split('/')[-1].split('.')[0]
    with open(fa) as f, open(sp_gene + '_set.fa', 'w') as outF:
        ID_seq_dict = {}
        seq_list = []
        for record in SeqIO.parse(f, 'fasta'):
            ID = record.id
            seq = str(record.seq)
            seq_list.append(seq)
            ID_seq_dict[seq] = ID
        for seqs in set(seq_list):
            IDs = ID_seq_dict[seqs]
            outF.write('>' + IDs + '\n' + seqs + '\n')


for fasta in glob.glob("*_set.fa"):
    sp = fasta.split('_')[0]
    sp_gene = fasta.split('_set')[0]
    f=open(sp_gene + "_set.fa","r");
    opened = False
    for line in f :
        if(line[0] == ">") :
            if(opened) :
                of.close()
            opened = True
            of=open(sp_gene + "_%s.fa" % (line[1:].rstrip()), "w")
            print(line[1:].rstrip())
        of.write(line)
    of.close()


sp_gene_list = []
for fasta in glob.glob("*gene*"):
    sp = fasta.split('_')[0]
    sp_gene_list.append(fasta)
    
def translate(sp_gene):
    spName = sp_gene.split('.')[0]
    print(sp_gene)
    unix('sixpack -sequence ' + sp_gene + ' -outfile ' + spName + '.sixpack -outseq ' + spName + '.sixpack.fa', shell=True)
    
for gene in sp_gene_list:
    translate(gene)


for AA_file in glob.glob("*.sixpack.fa"):
    peptide_list = []
    sp_gene = AA_file.split('_gene')[0]
    geneNum = AA_file.split('_gene')[-1].split('.')[0]
    sp = sp_gene.split('_')[0]
    gene = sp_gene.split('_')[1]
    chemo_ref_seq = chemo_gene_seq_dict[gene]
    with open(AA_file) as f, open(sp_gene + '.fasta', 'a+') as outF:
        seq_lens = []
        seq_len_dict = {}
        for record in SeqIO.parse(f, 'fasta'):
            seq = str(record.seq)
            seq_lens.append(len(seq))
            pep_len = str(len(seq))
            seq_len_dict[seq] = pep_len
            
        sp_gene_seq_idents = {}
        for pep, length in seq_len_dict.items():
            if int(length) >=100:
                alignments = aligner.align(chemo_ref_seq, pep)
                alignment = alignments[0]
                align_score = alignment.score
                sp_gene_seq_idents[pep] = align_score
        try:
            best_sp_gene_seq = max(sp_gene_seq_idents.values())
            for gene_seq, score in sp_gene_seq_idents.items():
                if score == best_sp_gene_seq:
                    outF.write('>' + sp + '|' + gene + geneNum + '\n' + gene_seq + '\n')
                    break
        except:
            continue
