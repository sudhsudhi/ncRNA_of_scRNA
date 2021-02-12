
**Instructions to run the python script:**

Make sure [pysam](https://pysam.readthedocs.io/en/latest/installation.html), [pandas](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html) and [numpy](https://pypi.org/project/numpy/) are installed.

To run the script, type the following command in your commandline (your current directory should be the one where the python script is saved):

```
python nocs_ncrna.py --bam_file path/to/bam_file/bam_file.bam --bai_file=path/to/bai_file/bai_file.bai --counts_csv path/to/csv_file/csv_file.csv --counts_csv_nocs  path/to/csv_file/csv_nocs_file.csv --pdf_file path/to/pdf_file/pdf_file.pdf --bed_file path/to/pdf_file/bed_file.bed

```

*Input files:*
Replace 'path/to/bam_file/bam_file.bam' and 'path/to/bai_file/bai_file.bai' with the paths to where the bam and bai files that you want to analyse are stored.

*Output files:*
The program generates four output files(a bed file of ncos-ncRNA (NOvel Cell-type Specific ncRNA) location in the genome, two csv files which stores the counts matrixes(one of genes and one of non-genes) and a pdf file that shows analytics). Replace 'path/to/csv_file/csv_file.csv', 'path/to/csv_file/csv_nocs_file.csv', 'path/to/pdf_file/pdf_file.pdf', 'path/to/pdf_file/bed_file.bed' with the paths to where you want these files to be generated at.(If in the current directory, just mention 'csv_file_name.csv')

If any of the 6 paths are not mentioned the program throws an error.

The command line option and the path to the sting can be separated either by a space or by "="( you could call the program as bam_file=/path/.. or bam_file /path/..)

You could also run using shorter commands as:
```
python nocs_ncrna.py -b path/to/bam_file/bam_file.bam -i path/to/bai_file/bai_file.bai -c path/to/csv_file/csv_file.csv -n path/to/csv_file/csv_nocs_file.csv -p path/to/pdf_file/pdf_file.pdf -d path/to/pdf_file/bed_file.bed

```
