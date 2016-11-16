### Analysis of public T-cell receptor repertoire

Datasets are available here:

* ``D1``: 500+ healthy donors from ImmunoSEQ (CDR3 beta amino acid sequences only): https://files.pub.cdr3.net/mikesh/SZtchNDr5FBEKXM_A6XjX3xU71mkJw-u/HIP.tar
* ``D2``: Our Rep-Seq data for various donors: https://files.pub.cdr3.net/mikesh/7xaduJSLM-aa4qJavkDTvSSOksxw0qMc/MILAB.tar **Unpublished data**

> NOTE: the same CDR3 amino acid sequence can be found multiple times in a given sample (due to convergent recombination), these should be de-duped.

### Tasks

1. Construct a database of public CDR3 beta amino acid sequence (clonotype) variants.
   1. Build the distribution of clonotype sequence occurrences in ``D1``. Test if the fraction of clonotypes found in X samples is distributed according to Zipf law. Look at super-public sequence found in majority of samples. Variants like ``CAVF`` are likely to be artefacts.
   2. Select variants found in at least 1%, 5% and 10% of samples (hereafter public clonotypes). For each public clonotype compute the fraction of samples in ``D2`` that contain it. Also compute the fraction of reads (i.e. sum of ``freq`` column values) for each public clonotype. Check the consistency (correlation) between public clonotype occurrence in ``D1`` and ``D2``. Select an optimal frequency threshold for defining a public clonotype from ``D1`` based on maximizing the correlation and number of occurrences in ``D2``.
2. Build co-occurrence graph for public clonotypes.
   1. Public clonotypes (nodes) should be joined with an edge if they are found together in at least 1 ``D1`` sample. 
   2. Edge weights are computed as the number of ``D1`` samples in which clonotypes are found together.
   3. Analyze the network properties of the graph. See [this](https://biodatamining.biomedcentral.com/track/pdf/10.1186/1756-0381-4-10?site=biodatamining.biomedcentral.com) paper for ideas.
3. Perform graph clustering based on methods described in [here](http://www.leonidzhukov.net/hse/2015/networks/papers/GraphClustering_Schaeffer07.pdf) ([implementation](http://micans.org/mcl/))
    1. Check if clusters have restricted Variable segment usage, i.e. a dominant Variable segment per cluster **[not available yet, provided upon ``D1`` update]**
    2. Check if clusters have distinct spectratype (distribution of CDR3 lengths).
    3. Check if CDR3 sequences within the same cluster are more similar than between clusters, e.g. compare average hamming distance
4. **[Data provided later]** Map donor cohorts with different HLA alleles onto public clonotype graph.