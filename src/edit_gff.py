#!/usr/bin/python3
import sys
import glob
from subprocess import call as unix


'''
Script to parse exonerate output to create
fasta files of coding sequences for all 
query genes in new species assembly
'''


#For each gff file, create temp file without unnecessary lines, run gff convert perl script, create fasta files of coding sequences using gffread
for gff in glob.glob("../*gff.filter"):
    sp = gff.split('/')[-1].split('_')[0]
    sp_gene = gff.split('/')[-1].split('.')[0]
    
    sp_gene_blast_hit = '../../blastout/sp_blast_extract/' + sp_gene + '.fasta'
    print(sp_gene)
    
    with open(gff) as f, open(gff + '.tmp', 'w') as outF:
        for line in f:
            if line.startswith('#'):
                continue
            elif line == '\n':
                continue
            else:
                lines = line.split('\t')
                if lines[2] != 'cds':
                    outF.write(line)
                    
                    
    unix('perl convert_exonerate_gff_to_std_gff.pl ' + gff + '.tmp '  + sp_gene + '.std.gff', shell=True)#Convert exonerate gff to standard gff files 

    unix('scp ../../blastout/sp_blast_extract/' + sp_gene + '.fasta .', shell=True)
    print('Creating fasta file...')
    unix('gffread -w ' + sp_gene + '.fa -g ' + sp_gene + '.fasta ' + sp_gene + '.std.gff', shell=True)#Run gffread to create fasta files
    
    print('Clearing dir..')
    unix('rm ' + sp_gene + '.fasta', shell=True)
    unix('rm *fai', shell=True)
    unix('rm ../*.tmp', shell=True)
    
    print('\n')
    
