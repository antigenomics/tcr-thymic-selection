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
As a TCR sequnces before selection we took the sequnces generated by OLGA $^1$ software. For naive cells after selection we took KECK dataset $^2$ - adaptive biotechnologies data (10 samples, approx. $10^6$ cdr3). Has read count of one = naive.
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
Kidera factors - the key fetures of peptides, first introduced by Kidera et al. $^3$. We analyzed each of 10 kidera factors genepairwise. Kideras with maximum FC were represented in Volcano plot above.
The most selected Kideras:  KF2, KF6 - Size; KF4 - Hydrophilicity; KF8 - Occurrence in alpha region
![newplot (28)](https://github.com/antigenomics/tcr-thymic-selection/assets/75636485/f1a9337a-0a24-48c7-9900-6354080907e6)

5. Genes clusters.  
We analysed gene clusters enriched in OLGA and rare in KECK and vise versa with VDJtools $^4$. Volcano plots for this clusters are presented below.
![image](https://github.com/antigenomics/tcr-thymic-selection/assets/75636485/2b8e9a05-b9ed-48e8-be27-63f67f7c2dc5)  

For all plots please see https://docs.google.com/presentation/d/12NM-7CLGjYuhLo4dROhbbfClH7IH4oHbi5svjMNscro/edit?usp=sharing

## Conclusions

1. C - strong negative selection
2. Nx[S,T] - glycosylation sites - strong negative selection
3. Other post-translational modification do not affect selection
4. Decrease in charge and increase in hydrophobicity. Extreme length are also negatively selected
5. Decrease in Kidera factors responsible for size and hydrophilicity. Increase in Kideras responsible for occurrence in alpha region
6. Negative selection toward “bad” aa
7. Clusters enriched in OLGA are short and R-reach
8. Clusters enriched in KECK are of normal size and G-reach
9. Clusters with glycosylation sites which survive after selection have J1-6 gene which is rare

## Methods
1. Generation of TCRs beta chain were carried out by OLGA software (1.2.4)
`olga-generate_sequences --humanTRB -n 10000000`
2. All statistics and visualisation except clusters enrichment were done in Python 3.7
3. Kidera factors and other phys-chemical properties were calculated in Peptides $^5$ (0.3.2)
4. Clusters enrichment comparison were carried out in VDJtools (1.2.1)  
`vdjtools --CalcDegreeStats` with default parameters

## References

1. Zachary Sethna et al., OLGA: fast computation of generation probabilities of B- and T-cell receptor amino acid sequences and motifs, Bioinformatics, Volume 35, Issue 17, 1 September 2019, Pages 2974–2981, https://doi.org/10.1093/bioinformatics/btz035
2. Emerson, R., DeWitt, W., Vignali, M. et al., Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire. Nat Genet 49, 659–665 (2017). https://doi.org/10.1038/ng.3822
3. Kidera, A., Konishi, Y., Oka, M. et al. Statistical analysis of the physical properties of the 20 naturally occurring amino acids. J Protein Chem 4, 23–55 (1985). https://doi.org/10.1007/BF01025492
4. Shugay M et al. VDJtools: Unifying Post-analysis of T Cell Receptor Repertoires. PLoS Comp Biol 2015; 11(11):e1004503-e1004503
5. Peptides lib: https://pypi.org/project/peptides/
