# -*- coding: utf-8 -*-
"""Processing for Genes PBMC.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1fzXypcVmoNMer6p-cYVmB02-UutgXzqC
"""

import pysam
import pandas as pd
from tqdm.auto import tqdm
import numpy as np
import itertools
import sys
import getopt
print("done importing")

#"--bam_file","--bai_file","--counts_csv","--pdf_file" are all necessarily rquired arguments.
# or you could pass as "-b","-i","-c","-p" respectively for the arguments
'''
sample execution of program:
python processing_for_genes.py -b /content/drive/My_Drive/pbmc_1k_protein_v3_possorted_genome_bam.bam -i /content/drive/My_Drive/pbmc_1k_protein_v3_possorted_genome_bam.bam.bai -c /content/drive/My_Drive/countss1.csv -p /content/drive/My_Drive/resultss1.pdf ' 

'''
def cmd_arguments()
  try:
        opts, args = getopt.getopt(sys.argv[1:],"b:i:c:p:",["bam_file=","bai_file=","counts_csv=","pdf_file="])
  except getopt.GetoptError:
        print('wrong argument format')
        print('run as : python processing_for_genes.py --bam_file path/to/bam_file/bam_file.bam --bai_file=path/to/bai_file/bai_file.bai --counts_csv path/to/csv_file/csv_file.csv --pdf_file path/to/pdf_file/pdf_file.pdf ')
        print('or run as : python processing_for_genes.py -b path/to/bam_file/bam_file.bam -i path/to/bai_file/bai_file.bai -c path/to/csv_file/csv_file.csv -p path/to/pdf_file/pdf_file.pdf ')
        
        print("sample execution of program:")
        print("python processing_for_genes.py -b /content/drive/My Drive/pbmc_1k_protein_v3_possorted_genome_bam.bam -i /content/drive/My Drive/pbmc_1k_protein_v3_possorted_genome_bam.bam.bai -c /content/drive/My Drive/countss1.csv -p /content/drive/My Drive/resultss1.pdf")
        sys.exit(2)

  '''code to confirm files exist and directories of csv and pdf exits'''

  for opt, arg in opts:
      if opt in ('-b', '--bam_file'):
          bam_file = arg
      elif opt in ('-i', '--bai_file'):
          bai_file = arg
      elif opt in ('-c', '--counts_csv'):
          counts_csv = arg
      elif opt in ('-p', '--pdf_file'):
          pdf_file = arg

  return bam_file,bai_file,counts_csv,pdf_file

def reader(bam_file,bai_file):
  samf = pysam.Samfile(bam_file, "rb")
  replicon_dict = dict([[replicon, {'seq_start_pos': 0,'seq_end_pos': length}] for replicon, length in zip(samf.references, samf.lengths)])
  #replicon_dict={'ref1': {'seq_start_pos': 0, 'seq_end_pos': 1}, 'ref2': {'seq_start_pos': 0, 'seq_end_pos': 2}, 'ref3': {'seq_start_pos': 0, 'seq_end_pos': 3}}

  samfile = pysam.AlignmentFile(bam_file, "rb", index_filename = bai_file)

  x=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
  list_tags = []

  for i in tqdm(range(0,len(x))):
      for read in samfile.fetch(x[i]):
          #x[i] is the value of the contig; contig in the sense of SAM files is the chromosome.
          #samfile.fetch(x[i]) gives an iterator over reads in x[i] chromosome.
          try:
              if read.has_tag("GX") and read.get_tag('NH')==1:
                  '''
                  GX= Semicolon-separated list of gene IDs that are compatible with this alignment. Gene IDs are specified with the gene_id key in the reference GTF attribute column.
                  NH= Number of reported alignments for query
                  CB=Cell barcode that is error-corrected and confirmed against a list of known-good barcode sequences
                  UB=UMI that is error-corrected among other molecular barcodes with the same cellular barcode and gene alignment
                  '''
                  list_tags.append([read.get_tag("CB"),read.get_tag("UB"),read.get_tag("GX")])
          except KeyError:
              continue


  list_tags_rm_dup_final = list(list_tags for list_tags,_ in itertools.groupby(list_tags))
  #list_tags_rm_dup_final=[[cb,ub,geneID1;geneID2;],[cb,ub,geneID1;geneID2;],...]

  return list_tags_rm_dup_final


  # '''
# # [k for k, g in groupby('AAAABBBCCDAABBB')] --> A B C D A B
# # [list(g) for k, g in groupby('AAAABBBCCD')] --> AAAA BBB CC D
# # recheck: exactly similar reads i.e same ub,cb,gx, will be reduced to one in this step; duplicates removed
# '''


bam_file,bai_file,counts_csv,pdf_file=cmd_arguments()

list_tags_rm_dup_final=reader(bam_file,bai_file)

print("done reading bam and bai files")

mul_genes=0
for i in list_tags_rm_dup_final:
  if ";" in i[2] or ":" in i[2]:
    mul_genes+=1
print("number of reads matching to multiple genes: "+str(mul_genes))

n=len(list_tags_rm_dup_final)


df_all_p1 = pd.DataFrame(list_tags_rm_dup_final[:int(n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p2 = pd.DataFrame(list_tags_rm_dup_final[int(n/10):int(2*n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p3 = pd.DataFrame(list_tags_rm_dup_final[int(2*n/10):int(3*n/10)], columns=["celltag","moltag","pseudoname"])
'''
df_all_p4 = pd.DataFrame(list_tags_rm_dup_final[int(3*n/10):int(4*n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p5 = pd.DataFrame(list_tags_rm_dup_final[int(4*n/10):int(5*n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p6 = pd.DataFrame(list_tags_rm_dup_final[int(5*n/10):int(6*n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p7 = pd.DataFrame(list_tags_rm_dup_final[int(6*n/10):int(7*n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p8 = pd.DataFrame(list_tags_rm_dup_final[int(7*n/10):int(8*n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p9 = pd.DataFrame(list_tags_rm_dup_final[int(8*n/10):int(9*n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p10 = pd.DataFrame(list_tags_rm_dup_final[int(9*n/10):], columns=["celltag","moltag","pseudoname"])
'''
print("done generating dataframes")


c1 = df_all_p1['celltag'].value_counts()
c2 = df_all_p2['celltag'].value_counts()
c3 = df_all_p3['celltag'].value_counts()
'''
c4 = df_all_p4['celltag'].value_counts()
c5 = df_all_p5['celltag'].value_counts()
c6 = df_all_p6['celltag'].value_counts()
c7 = df_all_p7['celltag'].value_counts()
c8 = df_all_p8['celltag'].value_counts()
c9 = df_all_p9['celltag'].value_counts()
c10 = df_all_p10['celltag'].value_counts()
'''
c11 = df_all_p1['pseudoname'].value_counts()
c22 = df_all_p2['pseudoname'].value_counts()
c33 = df_all_p3['pseudoname'].value_counts()
'''
c44 = df_all_p4['pseudoname'].value_counts()
c55 = df_all_p5['pseudoname'].value_counts()
c66 = df_all_p6['pseudoname'].value_counts()
c77 = df_all_p7['pseudoname'].value_counts()
c88 = df_all_p8['pseudoname'].value_counts()
c99 = df_all_p9['pseudoname'].value_counts()
c1010 = df_all_p10['pseudoname'].value_counts()
'''



df_all_p1_subset = df_all_p1[df_all_p1["celltag"].isin(c1[c1>0].index)]
#df_all_p1_subset is now a subset of df_all_p1 with a value in cell_tag column(read indexing with isin pandas)
df_all_p1_subset = df_all_p1_subset[df_all_p1_subset["pseudoname"].isin(c11[c11>0].index)]
#df_all_p1_subset is now a subset of df_all_p1 with a value in cell_tag colum AND a value in pseudoname column; i.e each entry has a cell barcode and is matched to a gene
df_all_p2_subset = df_all_p2[df_all_p2["celltag"].isin(c2[c2>0].index)]
df_all_p2_subset = df_all_p2_subset[df_all_p2_subset["pseudoname"].isin(c22[c22>0].index)]
df_all_p3_subset = df_all_p3[df_all_p3["celltag"].isin(c3[c3>0].index)]
df_all_p3_subset = df_all_p3_subset[df_all_p3_subset["pseudoname"].isin(c33[c33>0].index)]
'''
df_all_p4_subset = df_all_p4[df_all_p4["celltag"].isin(c4[c4>0].index)]
df_all_p4_subset = df_all_p4_subset[df_all_p4_subset["pseudoname"].isin(c44[c44>0].index)]
df_all_p5_subset = df_all_p5[df_all_p5["celltag"].isin(c5[c5>0].index)]
df_all_p5_subset = df_all_p5_subset[df_all_p5_subset["pseudoname"].isin(c55[c55>0].index)]
df_all_p6_subset = df_all_p6[df_all_p6["celltag"].isin(c6[c6>0].index)]
df_all_p6_subset = df_all_p6_subset[df_all_p6_subset["pseudoname"].isin(c66[c66>0].index)]
df_all_p7_subset = df_all_p7[df_all_p7["celltag"].isin(c7[c7>0].index)]
df_all_p7_subset = df_all_p7_subset[df_all_p7_subset["pseudoname"].isin(c77[c77>0].index)]
df_all_p8_subset = df_all_p8[df_all_p8["celltag"].isin(c8[c8>0].index)]
df_all_p8_subset = df_all_p8_subset[df_all_p8_subset["pseudoname"].isin(c88[c88>0].index)]
df_all_p9_subset = df_all_p9[df_all_p9["celltag"].isin(c9[c9>0].index)]
df_all_p9_subset = df_all_p9_subset[df_all_p9_subset["pseudoname"].isin(c99[c99>0].index)]
df_all_p10_subset = df_all_p10[df_all_p10["celltag"].isin(c10[c10>0].index)]
df_all_p10_subset = df_all_p10_subset[df_all_p10_subset["pseudoname"].isin(c1010[c1010>0].index)]
'''


counts_p1=df_all_p1_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p2=df_all_p2_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p3=df_all_p3_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
# .groupby(['pseudoname','celltag']).size() groups the rows with the same celltag AND matching to same pseudoname into a single group and stores the size of the group
# unstack converts a group into a datapoint by making the celltag a row_index and pseudoname a column index
# we get counts matrix at end of this
print(df_all_p1_subset.groupby(['pseudoname','celltag']).size())
print(counts_p1)  
'''
counts_p4=df_all_p4_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p5=df_all_p5_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p6=df_all_p6_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p7=df_all_p7_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p8=df_all_p8_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p9=df_all_p9_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p10=df_all_p10_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
'''


cell_bcode=set(counts_p1.columns).intersection(set(counts_p2.columns)).intersection(set(counts_p3.columns))
#.intersection(set(counts_p4.columns)).intersection(set(counts_p5.columns)).intersection(set(counts_p6.columns)).intersection(set(counts_p7.columns)).intersection(set(counts_p8.columns)).intersection(set(counts_p9.columns)).intersection(set(counts_p10.columns))
#I think this should be union
cell_bcode_2=set(counts_p1.columns).union(set(counts_p2.columns)).union(set(counts_p3.columns))

counts_p1 = counts_p1[cell_bcode]
counts_p2 = counts_p2[cell_bcode]
counts_p3 = counts_p3[cell_bcode]
'''
counts_p4 = counts_p4[cell_bcode]
counts_p5 = counts_p5[cell_bcode]
counts_p6 = counts_p6[cell_bcode]
counts_p7 = counts_p7[cell_bcode]
counts_p8 = counts_p8[cell_bcode]
counts_p9 = counts_p9[cell_bcode]
counts_p10 = counts_p10[cell_bcode]
'''
counts_full = pd.concat([counts_p1,counts_p2, counts_p3])
                         #counts_p4, counts_p5, counts_p6, counts_p7, counts_p8, counts_p9, counts_p10])
counts_new=counts_full.groupby(level=0, axis=0).sum()
print("done generating counts matrix")
#"/content/drive/My Drive/counts_genes_pbmc_10k_rmmulti_final.csv"
counts_new.to_csv(counts_csv)
print("done generating csv file")
num_genes=counts_new.astype(bool).sum(axis=1)
#num_genes is the number of genes that each cell matches to
num_genes
hist=num_genes.hist(bins=20)
for ax in hist.flatten():
    ax.set_xlabel("Number of genes each cell aligns to")
    ax.set_ylabel("Number of cells")


fig = hist.get_figure()
save_file_path="/content/drive/My Drive/"


fig.savefig(pdf_file)

