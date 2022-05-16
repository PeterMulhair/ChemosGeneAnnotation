import glob
from collections import Counter, defaultdict
from Bio import SeqIO

sp_genes = defaultdict(dict)
for fa in glob.glob("*.fasta"):
    sp_list = []
    gene = fa.split('_')[1].split('.')[0]
    sp = fa.split('_')[0]
    gene_count = 0
    gene_counts = {}
    with open(fa) as f:
        for record in SeqIO.parse(f, 'fasta'):
            gene_count+=1

    print(sp, gene, gene_count)
