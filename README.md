## Analysis of public T-cell receptor repertoire

The definition of public repertoire and public responses can be found in [Venturi et al. Nature Reviews Immunology 2008](http://www.nature.com.sci-hub.bz/nri/journal/v8/n3/pdf/nri2260.pdf)

The following large datasets were selected for studying public T-cell responses:

* ``D1``: 500+ healthy donors from ImmunoSEQ (CDR3 beta amino acid sequences only): ``D1/`` folder in repository root, this is a tar of gzipped lists of CDR3 amino acid sequences in each donor sample.
* ``D2``: Our Rep-Seq data for various donors: ``D2/`` folder in repository root contains datasets in VDJtools format. **This is an unpublished dataset!**. For more info on the file format check out [VDJtools documentation](http://vdjtools-doc.readthedocs.io/en/latest/input.html#input).

### Tasks

1. Construct a database of public CDR3 beta amino acid sequence (clonotype) variants.
   1. Build the distribution of clonotype sequence occurrences in ``D1``. Test if the fraction of clonotypes found in ``x`` samples is distributed according to [Zipf law](https://en.wikipedia.org/wiki/Zipf's_law). Look at super-public sequence found in majority of samples. Variants like ``CAVF`` are likely to be artefacts.
   2. Select variants found in at least 1%, 5% and 10% of samples (hereafter public clonotypes). For each public clonotype compute the fraction of samples in ``D2`` that contain it. Also compute the fraction of reads (i.e. sum of ``freq`` column values) for each public clonotype. Check the consistency (correlation) between public clonotype occurrence in ``D1`` and ``D2``. Select an optimal frequency threshold for defining a public clonotype from ``D1`` based on maximizing the correlation and number of occurrences in ``D2``.
2. Build co-occurrence graph for public clonotypes.
   1. Public clonotypes (nodes) should be joined with an edge if they are found together in at least 1 ``D1`` sample. 
   2. Edge weights are computed as the number of ``D1`` samples in which clonotypes are found together.
   3. The other way to compute edge weights and to filter unreliable edges is to use Hypergeometric model. If two clonotypes occur together in ``n_{ab}`` sample, each of them separately is found in ``n_{a}`` and ``n_{b}`` samples and there are ``n`` samples in total, the P-value can be computed as ``Hypgeom(n_{ab},n_{a},n_{b},n)``.
   4. Analyze the network properties of the graph. See [this](https://biodatamining.biomedcentral.com/track/pdf/10.1186/1756-0381-4-10?site=biodatamining.biomedcentral.com) paper for ideas.
3. Perform graph clustering based on methods described in [here](http://www.leonidzhukov.net/hse/2015/networks/papers/GraphClustering_Schaeffer07.pdf) ([implementation](http://micans.org/mcl/))
   1. Check if clusters have restricted Variable segment usage, i.e. a dominant Variable segment per cluster (note that V/J designations are only available for ``D2``)
   2. Check if clusters have distinct [spectratype](http://vdjtools-doc.readthedocs.io/en/latest/basic.html#id13) (distribution of CDR3 lengths).
   3. Check if CDR3 sequences within the same cluster are more similar than between clusters, e.g. compare average hamming distance.
4. **[Data provided later]** Map donor cohorts with different HLA alleles onto public clonotype graph. Our collaborators have HLA-typed repertoire samples + some are already published by [ImmunoSEQ](https://www.ncbi.nlm.nih.gov/pubmed/28002888)


> NOTE: the same CDR3 amino acid sequence can be found multiple times in a given sample (due to convergent recombination), these should be de-duped.
