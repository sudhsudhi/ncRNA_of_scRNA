

import pysam
import pandas as pd
from tqdm.auto import tqdm
import numpy as np
import itertools
import sys
import getopt
import matplotlib as plt
import matplotlib.pyplot
print("done importing")

#"--bam_file","--bai_file","--counts_csv","--pdf_file" are all necessarily rquired arguments.
# or you could pass as "-b","-i","-c","-p" respectively for the arguments
'''
sample execution of program:
python processing_for_genes.py -b /content/drive/My_Drive/pbmc_1k_protein_v3_possorted_genome_bam.bam -i /content/drive/My_Drive/pbmc_1k_protein_v3_possorted_genome_bam.bam.bai -c /content/drive/My_Drive/countss1.csv -p /content/drive/My_Drive/resultss1.pdf ' 

'''

try:
      opts, args = getopt.getopt(sys.argv[1:],"b:i:c:n:p:d:",["bam_file=","bai_file=","counts_csv=","counts_csv_nocs=","pdf_file=","bed_file="])
except getopt.GetoptError:
      print('wrong argument format')
      print('run as : python nocs_ncrna.py --bam_file path/to/bam_file/bam_file.bam --bai_file=path/to/bai_file/bai_file.bai --counts_csv path/to/csv_file/csv_file.csv --counts_csv_nocs  path/to/csv_file/csv_nocs_file.csv --pdf_file path/to/pdf_file/pdf_file.pdf --bed_file path/to/pdf_file/bed_file.bed')
      print("")
      print('or run as : python nocs_ncrna.py -b path/to/bam_file/bam_file.bam -i path/to/bai_file/bai_file.bai -c path/to/csv_file/csv_file.csv -n path/to/csv_file/csv_nocs_file.csv -p path/to/pdf_file/pdf_file.pdf -d path/to/pdf_file/bed_file.bed')
      print("")
      print("sample execution of program:")
      print("python nocs_ncrna.py -b /content/drive/My Drive/pbmc_1k_protein_v3_possorted_genome_bam.bam -i /content/drive/My Drive/pbmc_1k_protein_v3_possorted_genome_bam.bam.bai -c /content/drive/My Drive/countss1.csv -n /content/drive/My Drive/countss2.csv -p /content/drive/My Drive/resultss1.pdf -d /content/drive/My Drive/resultss_bed.bed")
      sys.exit(2)

'''code to confirm files exist and directories of csv and pdf exits'''

for opt, arg in opts:
    if opt in ('-b', '--bam_file'):
        bam_file = arg
    elif opt in ('-i', '--bai_file'):
        bai_file = arg
    elif opt in ('-c', '--counts_csv'):
        counts_csv = arg
    elif opt in ('-n', '--counts_csv_nocs'):
        counts_csv_nocs = arg
    elif opt in ('-p', '--pdf_file'):
        pdf_file = arg
    elif opt in ('-d', '--bed_file'):
        bed_file = arg

print(bed_file)

#bam_file = "/content/drive/My Drive/pbmc_1k_protein_v3_possorted_genome_bam.bam"

#bai_file = "/content/drive/My Drive/pbmc_1k_protein_v3_possorted_genome_bam.bam.bai"

print("done reading bam and bai files")

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
                GX=	Semicolon-separated list of gene IDs that are compatible with this alignment. Gene IDs are specified with the gene_id key in the reference GTF attribute column.
                NH=	Number of reported alignments for query
                CB=Cell barcode that is error-corrected and confirmed against a list of known-good barcode sequences
                UB=UMI that is error-corrected among other molecular barcodes with the same cellular barcode and gene alignment
                '''
                list_tags.append([read.get_tag("CB"),read.get_tag("UB"),read.get_tag("GX")])
        except KeyError:
            continue


list_tags_rm_dup_final = list(list_tags for list_tags,_ in itertools.groupby(list_tags))
# '''
# # [k for k, g in groupby('AAAABBBCCDAABBB')] --> A B C D A B
# # [list(g) for k, g in groupby('AAAABBBCCD')] --> AAAA BBB CC D
# # recheck: exactly similar reads i.e same ub,cb,gx, will be reduced to one in this step; duplicates removed
# '''

#list_tags_rm_dup_final=[[cb,ub,geneID1;geneID2;],[cb,ub,geneID1;geneID2;],...]
mul_genes=0
for i in list_tags_rm_dup_final:
  if ";" in i[2] or ":" in i[2]:
    mul_genes+=1
print("number of reads matching to multiple genes: "+str(mul_genes))

n=len(list_tags_rm_dup_final)


df_all_p1 = pd.DataFrame(list_tags_rm_dup_final[:int(n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p2 = pd.DataFrame(list_tags_rm_dup_final[int(n/10):int(2*n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p3 = pd.DataFrame(list_tags_rm_dup_final[int(2*n/10):int(3*n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p4 = pd.DataFrame(list_tags_rm_dup_final[int(3*n/10):int(4*n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p5 = pd.DataFrame(list_tags_rm_dup_final[int(4*n/10):int(5*n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p6 = pd.DataFrame(list_tags_rm_dup_final[int(5*n/10):int(6*n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p7 = pd.DataFrame(list_tags_rm_dup_final[int(6*n/10):int(7*n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p8 = pd.DataFrame(list_tags_rm_dup_final[int(7*n/10):int(8*n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p9 = pd.DataFrame(list_tags_rm_dup_final[int(8*n/10):int(9*n/10)], columns=["celltag","moltag","pseudoname"])
df_all_p10 = pd.DataFrame(list_tags_rm_dup_final[int(9*n/10):], columns=["celltag","moltag","pseudoname"])
print("done generating dataframes")


c1 = df_all_p1['celltag'].value_counts()
c2 = df_all_p2['celltag'].value_counts()
c3 = df_all_p3['celltag'].value_counts()
c4 = df_all_p4['celltag'].value_counts()
c5 = df_all_p5['celltag'].value_counts()
c6 = df_all_p6['celltag'].value_counts()
c7 = df_all_p7['celltag'].value_counts()
c8 = df_all_p8['celltag'].value_counts()
c9 = df_all_p9['celltag'].value_counts()
c10 = df_all_p10['celltag'].value_counts()
c11 = df_all_p1['pseudoname'].value_counts()
c22 = df_all_p2['pseudoname'].value_counts()
c33 = df_all_p3['pseudoname'].value_counts()
c44 = df_all_p4['pseudoname'].value_counts()
c55 = df_all_p5['pseudoname'].value_counts()
c66 = df_all_p6['pseudoname'].value_counts()
c77 = df_all_p7['pseudoname'].value_counts()
c88 = df_all_p8['pseudoname'].value_counts()
c99 = df_all_p9['pseudoname'].value_counts()
c1010 = df_all_p10['pseudoname'].value_counts()



df_all_p1_subset = df_all_p1[df_all_p1["celltag"].isin(c1[c1>0].index)]
#df_all_p1_subset is now a subset of df_all_p1 with a value in cell_tag column(read indexing with isin pandas)
df_all_p1_subset = df_all_p1_subset[df_all_p1_subset["pseudoname"].isin(c11[c11>0].index)]
#df_all_p1_subset is now a subset of df_all_p1 with a value in cell_tag colum AND a value in pseudoname column; i.e each entry has a cell barcode and is matched to a gene
df_all_p2_subset = df_all_p2[df_all_p2["celltag"].isin(c2[c2>0].index)]
df_all_p2_subset = df_all_p2_subset[df_all_p2_subset["pseudoname"].isin(c22[c22>0].index)]
df_all_p3_subset = df_all_p3[df_all_p3["celltag"].isin(c3[c3>0].index)]
df_all_p3_subset = df_all_p3_subset[df_all_p3_subset["pseudoname"].isin(c33[c33>0].index)]
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



counts_p1=df_all_p1_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p2=df_all_p2_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p3=df_all_p3_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
# .groupby(['pseudoname','celltag']).size() groups the rows with the same celltag AND matching to same pseudoname into a single group and stores the size of the group
# unstack converts a group into a datapoint by making the celltag a row_index and pseudoname a column index
# we get counts matrix at end of this

counts_p4=df_all_p4_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p5=df_all_p5_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p6=df_all_p6_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p7=df_all_p7_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p8=df_all_p8_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p9=df_all_p9_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p10=df_all_p10_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)



cell_bcode=set(counts_p1.columns).intersection(set(counts_p2.columns)).intersection(set(counts_p3.columns)).intersection(set(counts_p4.columns)).intersection(set(counts_p5.columns)).intersection(set(counts_p6.columns)).intersection(set(counts_p7.columns)).intersection(set(counts_p8.columns)).intersection(set(counts_p9.columns)).intersection(set(counts_p10.columns))
#I think this should be union
cell_bcode_2=set(counts_p1.columns).union(set(counts_p2.columns)).union(set(counts_p3.columns))

counts_p1 = counts_p1[cell_bcode]
counts_p2 = counts_p2[cell_bcode]
counts_p3 = counts_p3[cell_bcode]
counts_p4 = counts_p4[cell_bcode]
counts_p5 = counts_p5[cell_bcode]
counts_p6 = counts_p6[cell_bcode]
counts_p7 = counts_p7[cell_bcode]
counts_p8 = counts_p8[cell_bcode]
counts_p9 = counts_p9[cell_bcode]
counts_p10 = counts_p10[cell_bcode]

counts_full = pd.concat([counts_p1,counts_p2, counts_p3])
                         #counts_p4, counts_p5, counts_p6, counts_p7, counts_p8, counts_p9, counts_p10])
counts_new=counts_full.groupby(level=0, axis=0).sum()
print("done generating counts matrix")
#"/content/drive/My Drive/counts_genes_pbmc_10k_rmmulti_final.csv"
counts_new.to_csv(counts_csv)
print("done generating csv file")

num_genes=counts_new.astype(bool).sum(axis=1)
#num_genes is the number of genes that each cell matches to
print(num_genes)
fig1, axes = plt.pyplot.subplots(nrows=2, ncols=2)
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)
fig1.tight_layout()
hist=num_genes.hist(bins=20,ax=axes[0,0])
print("max_num_of _genes: "+str(num_genes.max()))
#fig = hist.get_figure()

#print(type(fig.axes[0]))
axes[0,0].set_xlabel("number of genes",fontsize=4)
axes[0,0].set_ylabel("number of cells",fontsize=4)
axes[0,0].set_title("Histogram of num of genes(top) and non-genes(bottom) that each cell matches to.",fontsize=6)

#fig.savefig(save_file_path+'figure.pdf')
#plt.pyplot.subplots_adjust(left=0.08)

num_genes.sort_values()
num_genes.plot(ax=axes[0,1])
axes[0,1].set_ylabel("number of genes",fontsize=7)
axes[0,1].set_xlabel("cell IDs",fontsize=7)
axes[0,1].set_xticks([])
axes[0,1].set_title("Num of genes that each cell matches to.",fontsize=6)

#------------------------------------
#Processing for non-genes:

#analysis for non_genes:
#read all non-genes
x=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
list_tags_all = []

for i in tqdm(range(0,len(x))):
    for read in samfile.fetch(x[i]):
        try:
            if (read.has_tag("GX")==False and read.get_tag("NH")==1):
                list_tags_all.append([read.get_tag("CB"),read.get_tag("UB"),"gene" + "-" + str(read.reference_name) + "-" +str(read.get_reference_positions()[0])+'-'+ str(read.get_reference_positions()[-1])])
        except KeyError:
            continue

list_tags_all_rm_dup = list(list_tags_all for list_tags_all,_ in itertools.groupby(list_tags_all))
#list_tags_all_rm_dup=[['GAGAAATAGAACCGCA-1', 'ACTGCTCTGAGT', 'gene-14-21509976-21510066'], ['GCACTAATCTTCTCAA-1', 'ACCTCCGTCGCA', 'gene-14-21510131-21510221']]

start=0
val=int(list_tags_all_rm_dup[start][2].split("-")[2])
collapse_len = 300

for i in tqdm(range(1,len(list_tags_all_rm_dup))):
    if (int(list_tags_all_rm_dup[i][2].split("-")[2]) - val>collapse_len):
        if (i-start>1):
            #update the inner elements
            s=list_tags_all_rm_dup[start][2].split("-")
            t=list_tags_all_rm_dup[i-1][2].split("-")
            new=s[0]+"-"+s[1]+"-"+s[2]+"-"+str(int(t[3]))
            for inner in list_tags_all_rm_dup[start:i]:
                inner[2] = new
        #update start to point to this pos
        start = i
        #update val to the val at this pos
        val = int(list_tags_all_rm_dup[i][2].split("-")[2])

list_tags_all_rm_dup_final = list(list_tags_all_rm_dup for list_tags_all_rm_dup,_ in itertools.groupby(list_tags_all_rm_dup))
n=len(list_tags_all_rm_dup_final)


df_all_p1 = pd.DataFrame(list_tags_all_rm_dup_final[:int(n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p2 = pd.DataFrame(list_tags_all_rm_dup_final[int(n/20):int(2*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p3 = pd.DataFrame(list_tags_all_rm_dup_final[int(2*n/20):int(3*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p4 = pd.DataFrame(list_tags_all_rm_dup_final[int(3*n/20):int(4*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p5 = pd.DataFrame(list_tags_all_rm_dup_final[int(4*n/20):int(5*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p6 = pd.DataFrame(list_tags_all_rm_dup_final[int(5*n/20):int(6*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p7 = pd.DataFrame(list_tags_all_rm_dup_final[int(6*n/20):int(7*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p8 = pd.DataFrame(list_tags_all_rm_dup_final[int(7*n/20):int(8*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p9 = pd.DataFrame(list_tags_all_rm_dup_final[int(8*n/20):int(9*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p10 = pd.DataFrame(list_tags_all_rm_dup_final[int(9*n/20):int(10*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p11 = pd.DataFrame(list_tags_all_rm_dup_final[int(10*n/20):int(11*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p12 = pd.DataFrame(list_tags_all_rm_dup_final[int(11*n/20):int(12*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p13 = pd.DataFrame(list_tags_all_rm_dup_final[int(12*n/20):int(13*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p14 = pd.DataFrame(list_tags_all_rm_dup_final[int(13*n/20):int(14*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p15 = pd.DataFrame(list_tags_all_rm_dup_final[int(14*n/20):int(15*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p16 = pd.DataFrame(list_tags_all_rm_dup_final[int(15*n/20):int(16*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p17 = pd.DataFrame(list_tags_all_rm_dup_final[int(16*n/20):int(17*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p18 = pd.DataFrame(list_tags_all_rm_dup_final[int(17*n/20):int(18*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p19 = pd.DataFrame(list_tags_all_rm_dup_final[int(18*n/20):int(19*n/20)], columns=["celltag","moltag","pseudoname"])
df_all_p20 = pd.DataFrame(list_tags_all_rm_dup_final[int(19*n/20):], columns=["celltag","moltag","pseudoname"])


c1 = df_all_p1['celltag'].value_counts()
c2 = df_all_p2['celltag'].value_counts()
c3 = df_all_p3['celltag'].value_counts()
c4 = df_all_p4['celltag'].value_counts()
c5 = df_all_p5['celltag'].value_counts()
c6 = df_all_p6['celltag'].value_counts()
c7 = df_all_p7['celltag'].value_counts()
c8 = df_all_p8['celltag'].value_counts()
c9 = df_all_p9['celltag'].value_counts()
c10 = df_all_p10['celltag'].value_counts()
c11 = df_all_p11['celltag'].value_counts()
c12 = df_all_p12['celltag'].value_counts()
c13 = df_all_p13['celltag'].value_counts()
c14 = df_all_p14['celltag'].value_counts()
c15 = df_all_p15['celltag'].value_counts()
c16 = df_all_p11['celltag'].value_counts()
c17 = df_all_p12['celltag'].value_counts()
c18 = df_all_p13['celltag'].value_counts()
c19 = df_all_p14['celltag'].value_counts()
c20 = df_all_p15['celltag'].value_counts()
c111 = df_all_p1['pseudoname'].value_counts()
c22 = df_all_p2['pseudoname'].value_counts()
c33 = df_all_p3['pseudoname'].value_counts()
c44 = df_all_p4['pseudoname'].value_counts()
c55 = df_all_p5['pseudoname'].value_counts()
c66 = df_all_p6['pseudoname'].value_counts()
c77 = df_all_p7['pseudoname'].value_counts()
c88 = df_all_p8['pseudoname'].value_counts()
c99 = df_all_p9['pseudoname'].value_counts()
c1010 = df_all_p10['pseudoname'].value_counts()
c1111 = df_all_p11['pseudoname'].value_counts()
c1212 = df_all_p12['pseudoname'].value_counts()
c1313 = df_all_p13['pseudoname'].value_counts()
c1414 = df_all_p14['pseudoname'].value_counts()
c1515 = df_all_p15['pseudoname'].value_counts()
c1616 = df_all_p16['pseudoname'].value_counts()
c1717 = df_all_p17['pseudoname'].value_counts()
c1818 = df_all_p18['pseudoname'].value_counts()
c1919 = df_all_p19['pseudoname'].value_counts()
c2020 = df_all_p20['pseudoname'].value_counts()


df_all_p1_subset = df_all_p1[df_all_p1["celltag"].isin(c1[c1>20].index)]
df_all_p1_subset = df_all_p1_subset[df_all_p1_subset["pseudoname"].isin(c111[c111>20].index)]
df_all_p2_subset = df_all_p2[df_all_p2["celltag"].isin(c2[c2>20].index)]
df_all_p2_subset = df_all_p2_subset[df_all_p2_subset["pseudoname"].isin(c22[c22>20].index)]
df_all_p3_subset = df_all_p3[df_all_p3["celltag"].isin(c3[c3>20].index)]
df_all_p3_subset = df_all_p3_subset[df_all_p3_subset["pseudoname"].isin(c33[c33>20].index)]
df_all_p4_subset = df_all_p4[df_all_p4["celltag"].isin(c4[c4>20].index)]
df_all_p4_subset = df_all_p4_subset[df_all_p4_subset["pseudoname"].isin(c44[c44>20].index)]
df_all_p5_subset = df_all_p5[df_all_p5["celltag"].isin(c5[c5>20].index)]
df_all_p5_subset = df_all_p5_subset[df_all_p5_subset["pseudoname"].isin(c55[c55>20].index)]
df_all_p6_subset = df_all_p6[df_all_p6["celltag"].isin(c6[c6>20].index)]
df_all_p6_subset = df_all_p6_subset[df_all_p6_subset["pseudoname"].isin(c66[c66>20].index)]
df_all_p7_subset = df_all_p7[df_all_p7["celltag"].isin(c7[c7>20].index)]
df_all_p7_subset = df_all_p7_subset[df_all_p7_subset["pseudoname"].isin(c77[c77>20].index)]
df_all_p8_subset = df_all_p8[df_all_p8["celltag"].isin(c8[c8>20].index)]
df_all_p8_subset = df_all_p8_subset[df_all_p8_subset["pseudoname"].isin(c88[c88>20].index)]
df_all_p9_subset = df_all_p9[df_all_p9["celltag"].isin(c9[c9>20].index)]
df_all_p9_subset = df_all_p9_subset[df_all_p9_subset["pseudoname"].isin(c99[c99>20].index)]
df_all_p10_subset = df_all_p10[df_all_p10["celltag"].isin(c10[c10>20].index)]
df_all_p10_subset = df_all_p10_subset[df_all_p10_subset["pseudoname"].isin(c1010[c1010>20].index)]
df_all_p11_subset = df_all_p11[df_all_p11["celltag"].isin(c11[c11>20].index)]
df_all_p11_subset = df_all_p11_subset[df_all_p11_subset["pseudoname"].isin(c1111[c1111>20].index)]
df_all_p12_subset = df_all_p12[df_all_p12["celltag"].isin(c12[c12>20].index)]
df_all_p12_subset = df_all_p12_subset[df_all_p12_subset["pseudoname"].isin(c1212[c1212>20].index)]
df_all_p13_subset = df_all_p13[df_all_p13["celltag"].isin(c13[c13>20].index)]
df_all_p13_subset = df_all_p13_subset[df_all_p13_subset["pseudoname"].isin(c1313[c1313>20].index)]
df_all_p14_subset = df_all_p14[df_all_p14["celltag"].isin(c14[c14>20].index)]
df_all_p14_subset = df_all_p14_subset[df_all_p14_subset["pseudoname"].isin(c1414[c1414>20].index)]
df_all_p15_subset = df_all_p15[df_all_p15["celltag"].isin(c15[c15>20].index)]
df_all_p15_subset = df_all_p15_subset[df_all_p15_subset["pseudoname"].isin(c1515[c1515>20].index)]
df_all_p16_subset = df_all_p16[df_all_p16["celltag"].isin(c16[c16>20].index)]
df_all_p16_subset = df_all_p16_subset[df_all_p16_subset["pseudoname"].isin(c1616[c1616>20].index)]
df_all_p17_subset = df_all_p17[df_all_p17["celltag"].isin(c17[c17>20].index)]
df_all_p17_subset = df_all_p17_subset[df_all_p17_subset["pseudoname"].isin(c1717[c1717>20].index)]
df_all_p18_subset = df_all_p18[df_all_p18["celltag"].isin(c18[c18>20].index)]
df_all_p18_subset = df_all_p18_subset[df_all_p18_subset["pseudoname"].isin(c1818[c1818>20].index)]
df_all_p19_subset = df_all_p19[df_all_p19["celltag"].isin(c19[c19>20].index)]
df_all_p19_subset = df_all_p19_subset[df_all_p19_subset["pseudoname"].isin(c1919[c1919>20].index)]
df_all_p20_subset = df_all_p20[df_all_p20["celltag"].isin(c20[c20>20].index)]
df_all_p20_subset = df_all_p20_subset[df_all_p20_subset["pseudoname"].isin(c2020[c2020>20].index)]

counts_p1=df_all_p1_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p2=df_all_p2_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p3=df_all_p3_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p4=df_all_p4_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p5=df_all_p5_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p6=df_all_p6_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p7=df_all_p7_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p8=df_all_p8_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p9=df_all_p9_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p10=df_all_p10_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p11=df_all_p11_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p12=df_all_p12_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p13=df_all_p13_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p14=df_all_p14_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p15=df_all_p15_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p16=df_all_p16_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p17=df_all_p17_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p18=df_all_p18_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p19=df_all_p19_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)
counts_p20=df_all_p20_subset.groupby(['pseudoname','celltag']).size().unstack('celltag', fill_value=0)

cell_bcode=set(counts_p1.columns).intersection(set(counts_p2.columns)).intersection(set(counts_p3.columns)).intersection(set(counts_p4.columns)).intersection(set(counts_p5.columns)).intersection(set(counts_p6.columns)).intersection(set(counts_p7.columns)).intersection(set(counts_p8.columns)).intersection(set(counts_p9.columns)).intersection(set(counts_p10.columns)).intersection(set(counts_p11.columns)).intersection(set(counts_p12.columns)).intersection(set(counts_p13.columns)).intersection(set(counts_p14.columns)).intersection(set(counts_p15.columns)).intersection(set(counts_p16.columns)).intersection(set(counts_p17.columns)).intersection(set(counts_p18.columns)).intersection(set(counts_p19.columns)).intersection(set(counts_p20.columns))

counts_p1 = counts_p1[cell_bcode]
counts_p2 = counts_p2[cell_bcode]
counts_p3 = counts_p3[cell_bcode]
counts_p4 = counts_p4[cell_bcode]
counts_p5 = counts_p5[cell_bcode]
counts_p6 = counts_p6[cell_bcode]
counts_p7 = counts_p7[cell_bcode]
counts_p8 = counts_p8[cell_bcode]
counts_p9 = counts_p9[cell_bcode]
counts_p10 = counts_p10[cell_bcode]
counts_p11 = counts_p11[cell_bcode]
counts_p12 = counts_p12[cell_bcode]
counts_p13 = counts_p13[cell_bcode]
counts_p14 = counts_p14[cell_bcode]
counts_p15 = counts_p15[cell_bcode]
counts_p16 = counts_p16[cell_bcode]
counts_p17 = counts_p17[cell_bcode]
counts_p18 = counts_p18[cell_bcode]
counts_p19 = counts_p19[cell_bcode]
counts_p20 = counts_p20[cell_bcode]

counts_full_ng = pd.concat([counts_p1,counts_p2, counts_p3, counts_p4, counts_p5, counts_p6, counts_p7, counts_p8, counts_p9, counts_p10, counts_p11,counts_p12, counts_p13, counts_p14, counts_p15, counts_p16, counts_p17, counts_p18, counts_p19, counts_p20])
counts_new_ng=counts_full_ng.groupby(level=0, axis=0).sum()
counts_new_ng.to_csv(counts_csv_nocs)




num_cells_ng=counts_new_ng.astype(bool).sum(axis=1)
#num_cells_ng is the number of cells that each gene matches to
num_count_ng=counts_new_ng.sum(axis=1)
#num_count_ng is the total number of reads that aligned with the TAR
median_expression=counts_new_ng.replace(0,np.nan).median(axis=1,skipna=True)
#median_expression is the median of the non-zero values in each row

print("num_cells_ng:")
print(num_cells_ng)
print("num_count_ng:")
print(num_count_ng)
print("median_expression:")
print(median_expression)

fin=pd.concat([num_cells_ng, num_count_ng,median_expression], axis=1)
fin.sort_values(by=0, inplace=True, ascending=False) #sorted based on number of cells that each TAR matches to
print("BED FILE:")
print(fin)
 
f=open(bed_file,"w")
#list_tags_all_rm_dup_final=[['GAGAAATAGAACCGCA-1', 'ACTGCTCTGAGT', 'gene-14-21509976-21510066'], ['GCACTAATCTTCTCAA-1', 'ACCTCCGTCGCA', 'gene-14-21510131-21510221']]
n=0
f.write("#chrom"+"\t"+ "chromStart"+"\t"+ "chromEnd"+"\t"+"num_cells_aligned_to"+"\t"+"num_molecules_aligned_to"+"\t"+"median_of_nonzero_count_vals"+"\n")
for index, row in fin.iterrows():
    #index=gene-14-21509976-21510066
    k=index.split("-")
    f.write("chr"+k[1]+"\t"+k[2]+"\t"+k[3]+"\t"+str(row[0])+"\t"+str(row[1])+"\t"+str(row[2])+"\n")
    if n<10:print(index,row[0],row[1],row[2])
    n+=1
f.close()

hist_ng=num_cells_ng.hist(bins=20,ax=axes[1,0],figure=fig1)
#print("max_num_of _genes: "+str(num_genes.max()))
#print(type(fig.axes[0]))
axes[1,0].set_xlabel("number of non genes",fontsize=7)
axes[1,0].set_ylabel("number of cells",fontsize=7)
#axes[1,0].set_title("Histogram of number of non-genes that each cell matches to.",fontsize=7)


num_cells_ng.sort_values()
num_cells_ng.plot(ax=axes[1,1],figure=fig1)
axes[1,1].set_ylabel("number of non genes",fontsize=7)
axes[1,1].set_xlabel("cell IDs",fontsize=7)
axes[1,1].set_xticks([])
axes[1,1].set_title("number of non-genes that each cell matches to.",fontsize=7)
fig1.savefig(pdf_file)