# False Duplication Identification

A script for false duplication identification with a simple example (a chromosome sequence in zebra finch).
 
NOTE: This script was not made for extensive using. This script have focused on the specific data set of assemblies. Many processed data (e.g. .maf, .bam, .meryl, purge_dups.bed..) with assemblies should be prepared before to run the script.

This script find the false duplications in the assemblies based on the alignment between two assemblies with depth-coverage, assembly gap, read-depth gap and paired read presence/absence validation.
The duplications found in self-alignment by purge_dups can be combined and merged to the result. 


Reference: Ko, Byung June, et al. "Widespread false gene gains caused by duplication errors in genome assemblies." bioRxiv (2021).

## Dependency
- meryl v1.0
- samtools
- bedtools
- pickle (python3)
- numpy (python3)
- parmap (python3)
- subprocess (python3)
- cigar (python3)
- multiprocessing (python3)

## The list of input files to prepare
### 1. Multiple alignment format (.maf) file
  - .maf file generated by Hal (https://github.com/ComparativeGenomicsToolkit/hal) from .hal file.
  - .hal generated by Cactus alignment (https://github.com/ComparativeGenomicsToolkit/cactus) of both assemblies.
### 2. Purge_dups.bed
  - .bed file generated by purge_dups package in the addloc branch of github (https://github.com/dfguan/purge_dups/tree/add_loc).
### 3. K-mer DB
  - K-mer DB generated by meryl v1.0 (v1.3 is not available in this script).
  - Both two assemblies and read DB have to be generated. Follow the process in Merqury (https://github.com/marbl/merqury/wiki/1.-Prepare-meryl-dbs).
### 4. Read mapping file
  - .bam for reads to be mapped to assembly (paired-end reads have to be mapped for this script).
### 5. Pairwise mapping format (.paf) file
  - The result of self-alignment in purge_dups process (https://github.com/dfguan/purge_dups).

### 6. Assembly gap file
  - Assembly gap position information in NCBI ftp.
  
### 7. Chromosome name conversion (.txt) file
  - Tab-deliminated scaffold name of assembly gap file and fasta (e.g. CM_111111 NC_111111) used for conversion of chromosome/scaffold name in the step of assembly gap parsing.
  - The scaffold name in first column changed to the name in second column.
  - If you don't need to conversion (i.e. if the scaffold names in the assembly gap file and the fasta file are same), set the second column as same as first column (e.g. CM_111111  CM_111111).
  - Put the name of both assemblies in a file. All scaffold name in alignment file (.maf) have to be included.
 

## Parameter setting in the configure file
- Cores: Set the processor used in the work. 
- HapChr: The part of genome of haploid (heterozygous sex chromosomes or mitochondiral DNA) may be specified as the name of scaffolds with tab-deliminated in 'HapChr'(e.g. HapChr Zeb_X1  Zeb_X2  Zeb_Y1 Zeb_Y2).
- InsertSizeCutoff: Set the insert size of the paired read (This information will be used for discordant read finding).
- Ref, Tar: The name of the assemblies in .maf.
- Asm*DepthCutoff: Depth coverage threshold to identify false duplication. The depth coverage under the AsmDepthCutoff is the potential false duplication.
- Asm*HapDepthCutoff: Depth coverage threshold to identify false duplication for haploid genome. The depth coverage under the AsmHapDepthCutoff is the potential false duplication.
- Asm*SECutoff: Depth coverage threshold to identify regions of assembly generated by sequencing error. 

## Run
```
cd FalseDuplication
python3 FDFinder.py <configure.txt>
```
<> : required

- The letters of first column in the configure file must be kept as same in the example file to role as the identifier of each parameter.
- Asm*FD_merged.bed will be generated if there was no errors during the run.
- Asm*FD.bed can be refered to capture the loci information of the sister sequence of false duplication.


## Example files
- Example.tar.gz contains example files (except k-mer DBs) to run the script. The files of zebra ficnh Chromosome 25 of previous assembly (GCF_000151805.1) and recent assembly (GCF_003957565.1) are included. 
- External links of K-mer DBs (two assembly DB, one read DB) for the example:  
https://drive.google.com/file/d/1bXGovqdAHQJbh2E1nozb_eeamB4iQNho/view?usp=sharing
https://drive.google.com/file/d/1jGo9O93etp-FdNYuOm1tXigh714kG0bY/view?usp=sharing
https://drive.google.com/file/d/1_qX3chSwVujlD-Jg0B6xgLvjHQnvXGzU/view?usp=sharing
- After decompressing the Example.tar.gz and k-merDBs, set the path of each file/folder in Configure.txt
