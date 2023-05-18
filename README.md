# Thymic selection proj

## Aims:
1. Study thymic selection using rep-seq data

## Objectives:
1. Generate TCR data using OLGA software and use it as TCR data before thymic selection
2. Compare gene usage between OLGA and experimental sample
3. Compare k-mers usage -//-
4. Study changes in chemo-physical properties of cdr3 overall and gene-pairwise
5. Study selected clusters of CDR3 sequences

## Results:

0. Data generation.
As a TCR sequnces before selection we took the sequnces generated by OLGA[] software. For naive cells after selection we took KECK dataset[] - adaptive biotechnologies data (10 samples, approx. 10^6 cdr3). Has read count of one = naive.
1. Gene usages plots.  
We calculated the genes usage frequences before and after thymic selection and compared it:

![image](https://github.com/antigenomics/tcr-thymic-selection/assets/75636485/ba8ab273-2952-45cc-bf16-51c7e9ba34c4)

![image](https://github.com/antigenomics/tcr-thymic-selection/assets/75636485/e9058d05-1c3f-4413-a633-f8a0b07983ab)

2. Kmers sequences. 
We compare the k-mers frequences between OLGA and KECK samples. Maximum attention were paid to 1-mers (single aa) and 3-mers.  
1-mers:  
![image](https://github.com/antigenomics/tcr-thymic-selection/assets/75636485/5cffcfa5-1ce3-430d-ab0d-31b300e802cb)  
3-mers:  
![image](https://github.com/antigenomics/tcr-thymic-selection/assets/75636485/6920d633-6b7f-4961-8813-4e02d5eb8744)

3. General chemo-physical proeprties.
We calculared key chemo-physical proeprties (length, charge, hydrophobicity) for all OLGA and KECK TCRs.  
![image](https://github.com/antigenomics/tcr-thymic-selection/assets/75636485/4d6fefb9-c5ae-43dd-aa64-941c7689b366)
![image](https://github.com/antigenomics/tcr-thymic-selection/assets/75636485/42be7cf5-b2c4-40fd-abdb-f6176dba4561)
![image](https://github.com/antigenomics/tcr-thymic-selection/assets/75636485/aeee9554-6244-4c86-8f06-75f0e407ec73)

4. Kidera factors analysis.
Kidera factors - the key fetures of peptides, first introduced by Kidera et al[]. We analyzed each of 10 kidera factors genepairwise. Kideras with maximum FC were represented in Volcano plot above.

![newplot (28)](https://github.com/antigenomics/tcr-thymic-selection/assets/75636485/f1a9337a-0a24-48c7-9900-6354080907e6)

5. Genes clusters.  
We analysed gene clusters enriched in OLGA and rare in KECK and vise versa. Volcano plots for this clusters are presented below.
![image](https://github.com/antigenomics/tcr-thymic-selection/assets/75636485/2b8e9a05-b9ed-48e8-be27-63f67f7c2dc5)
