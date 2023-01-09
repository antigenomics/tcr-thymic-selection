### Thymic selection proj

1. Generate pre-selection TCR repertoire from IGoR/OLGA software using model trained on Adaptive data
2. Get adaptive data for naive CD4/8 T-cell TCRs
3. Compare based on basic characteristics: Variable segment usage, CDR3 length, etc.
4. Check 8000 amino acid 3-mers based on over/under-representation (hypergeom+fold = volcano plots) in pre- and post-selection repertoire, see if some K-mers are positively or negatively selected.
5. Run ALICE/TCRNET analysis using a) pre-selection as control, CD4/CD8 as subject sample, b) CD4/CD8 as control and pre-selection as subject sample = negaive selection.
6. Extract subgraphs enriched in neighbours from 5. - see what motifs we get, are they similar/distinct between CD4/CD8, what are the properties (V,CDR3 len, hydroph/charged/.. AAs, etc) of TCRs in those motifs.
7. Check if we can find these motifs somewhere in our cancer/autoimmunity/antigen-specificity data.