#!/usr/bin/python3
import sys
import argparse
from collections import OrderedDict, defaultdict


'''
Exonerate output parser to get minimal set of 
annotated genes, filtering overlapping genes
and taking the longest ORF
'''


parse = argparse.ArgumentParser()

parse.add_argument("--gff",type=str, help="name of gff file to parse",required=True)

args = parse.parse_args()


sp = args.gff.split('_')[0]

#Create dictionary of all gene ORFs in gff output
gene_ranges_list = []
gene_range_remove = []
gene_exon_lens = {}
contig_gene_range = defaultdict(list)
with open(args.gff) as f:
    lines = f.read()
    hits = lines.split('# --- END OF GFF DUMP ---')
    for hit in hits:
        gene_intron_lens = defaultdict(list)
        intron_lens = []
        exon_lens = []
        exon_count = 0
        hit_info = hit.split('\n')
        for info in hit_info:
            if info.startswith('#'):
                continue
            else:
                gene_info = info.split('\t')
                if len(gene_info) > 1:
                    if gene_info[2] == 'intron':
                        intron_length = len(range(int(gene_info[3]), int(gene_info[4])))
                        intron_lens.append(intron_length)
                        

                    if gene_info[2] == 'gene':
                        if int(gene_info[3]) < int(gene_info[4]):
                            gene_start = gene_info[3]
                            gene_end = gene_info[4]
                        else:
                            gene_end = gene_info[3]
                            gene_start = gene_info[4]
                            
                        contig = gene_info[0]
                            
                        gene_range = range(int(gene_start), int(gene_end))
                        contig_gene_range[contig].append(gene_range)
                        gene_ranges_list.append(gene_range)


#For each contig with an annotated genes, filter overlapping genes retaining the largest ORF, and write to a new gff file        
total_gene_list = []
total_subset_gene_list = []
exonerate_hit_set = []
for contig, gene_ranges in contig_gene_range.items():
    
    #print(contig)
    
    gene_ranges_dict = {}
    for gene_range in set(gene_ranges):
        gene_ranges_dict[len(gene_range)] = gene_range

    ordered_range_list = []
    gene_ranges_dict = OrderedDict(sorted(gene_ranges_dict.items()))
    for k, v in reversed(gene_ranges_dict.items()):
        ordered_range_list.append(v)


    subset_gene_list = []
    filter_ordered_range_list = []
    for ranges in ordered_range_list:
        if ranges not in subset_gene_list:
            filter_ordered_range_list.append(ranges)

            
    for i in range(0, len(filter_ordered_range_list)):
        longest_gene_index = -1
        overlap = False
        range_i = filter_ordered_range_list[i]
        #    print(i, range_i)

        for j in range(0, len(filter_ordered_range_list)):
            range_j = filter_ordered_range_list[j]
            if range_i != range_j:        
                set_range_i = set(range_i)
                
                if set_range_i.intersection(range_j):
                    overlap = True
                    range_longest = max(range_i, range_j, key=len)
                    if range_longest not in subset_gene_list:
                        subset_overlap = False
                        for ranges in subset_gene_list:
                            if set(range_longest).intersection(ranges):
                                subset_overlap = True
                        if subset_overlap == False:
                    
                            subset_gene_list.append(range_longest)

        if overlap == False:
            subset_gene_list.append(range_i)
            
    #print(sorted(subset_gene_list, key=lambda r: r.start))
    #print(args.gff, len(subset_gene_list))

    total_gene_list.append(len(subset_gene_list))



    for sub_gene_ranges in sorted(subset_gene_list, key=lambda r: r.start):
        
        sub_gene_range_start = sub_gene_ranges[0]
        sub_gene_range_end = sub_gene_ranges[-1] + 1

        with open(args.gff) as f1:
            exoner_lines = f1.read()
            exoner_hits = exoner_lines.split('# --- END OF GFF DUMP ---')
            for exoner_hit in exoner_hits:
                exoner_hit_info = exoner_hit.split('\t' + 'similarity')
                try:
                    exoner_gene_info = exoner_hit_info[0].split('\t')
                    exoner_contig = exoner_gene_info[0].split('\n')[-1]
                    exoner_gene_start = exoner_gene_info[3].strip()
                    exoner_gene_end = exoner_gene_info[4].strip()
                except:
                    continue

                if contig == exoner_contig:
                    #print(sub_gene_range_start, sub_gene_range_end)
                    #print(exoner_gene_start, exoner_gene_end)
                    if (int(sub_gene_range_start) == int(exoner_gene_start)) and (int(sub_gene_range_end) == int(exoner_gene_end)):
                        exonerate_hit_set.append(exoner_hit)

                        break


print(args.gff, sum(total_gene_list))
#print('\n')

#with open('chemo_summary.tsv', 'a+') as outF:
#    outF.write(args.gff.split('.')[0] + '\t' + str( sum(total_gene_list)) + '\n')

outF = open(args.gff + '.filter', 'w')
for filter_hit in set(exonerate_hit_set):
    outF.write(filter_hit)
outF.close()


