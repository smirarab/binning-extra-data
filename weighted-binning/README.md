# Weighted Binning

Bayzid, Md. Shamsuzzoha, Siavash Mirarab, Bastien Boussau, and Tandy Warnow. “Weighted Statistical Binning: Enabling Statistically Consistent Genome-Scale Phylogenetic Analyses.” PLoS ONE 10, no. 6 (January 18, 2015): e0129183. doi:10.1371/journal.pone.0129183.

* Most of the datasets used in this study are available through [../main](../main) directory. 

* The new datasets generated for this study are available on [figshare, with DOI: 10.6084/m9.figshare.1411146](http://dx.doi.org/10.6084/m9.figshare.1411146). 

Below is the list of files on FigShare:

### [10-taxon datasets](https://figshare.com/articles/dataset/plos_datasets/1411146?file=3331325): 20 replicates, 200 genes, 100bp

File contents:

* Perl script for simulating indelible alignments: `post_stidsim.pl`
* True species trees (Higher ILS): `10-taxon/higher-ILS/true-speciestrees/R<X>.true.tre` - where X is the replicate number (1...20).
* True species trees (Lower ILS): `10-taxon/lower-ILS/true-speciestrees/R<X>.true.tre` - where X is the replicate number (1...20).
* True gene tree (Higer ILS): `10-taxon/higher-ILS/true-genetrees/R<X>/<Y>/true.gt` - where X is the replicate number (1...20), and Y is the gene id (0001...0200).
* True gene trees (Lower ILS): `10-taxon/lower-ILS/true-genetrees/R<X>/<Y>/true.gt` - where X is the replicate number (1...20), and Y is the gene id (0001...0200).
* Gene alignments (Higer ILS): `10-taxon/higher-ILS/estimated-genetrees/R<X>/<Y>/truegene.phy` - where X is the replicate number (1...20), and Y is the gene id (0001...0200).
* Gene alignments (Lower ILS): `10-taxon/lower-ILS/estimated-genetrees/R<X>/<Y>/truegene.phy` - where X is the replicate number (1...20), and Y is the gene id (0001...0200).
* Estimated BestML gene trees (Higher ILS): `10-taxon/higher-ILS/estimated-genetrees/R<X>/<Y>/RAxML_bipartitions.final.f200` - where X is the replicate number (1...20), and Y is the gene id (0001...0200).
* Estimated bootstrap gene trees (Higher ILS): `10-taxon/higher-ILS/estimated-genetrees/R<X>/<Y>/RAxML_bootstrap.all` - where X is the replicate number (1...20), and Y is the gene id (0001...0200).
* Estimated BestML gene trees (Lower ILS): `10-taxon/lower-ILS/estimated-genetrees/R<X>/<Y>/RAxML_bipartitions.final.f200` - where X is the replicate number (1...20), and Y is the gene id (0001...0200).
* Estimated bootstrap gene trees (Lower ILS): `10-taxon/lower-ILS/estimated-genetrees/R<X>/<Y>/RAxML_bootstrap.all` - where X is the replicate number (1...20), and Y is the gene id (0001...0200).
* Estimated BestML supergene trees (Higher ILS): `10-taxon/higher-ILS/estimated-supergenetrees/supergenes-<T>/R<X>/bin.<Z>.txt/RAxML_bipartitions.final.f200` - where X is the replicate number (1...20), T is binning threshold (50 or 75), and Z is the bin id.
* Estimated bootstrap supergene trees (Higher ILS): `10-taxon/higher-ILS/estimated-supergenetrees/supergenes-<T>/R<X>/bin.<Z>.txt/RAxML_bootstrap.all` - where X is the replicate number (1...20), T is binning threshold (50 or 75), and Z is the bin id.
* Estimated BestML supergene trees (Lower ILS): `10-taxon/lower-ILS/estimated-supergenetrees/supergenes-<T>/R<X>/bin.<Z>.txt/RAxML_bipartitions.final.f200` - where X is the replicate number (1...20), T is binning threshold (50 or 75), and Z is the bin id.
* Estimated bootstrap supergene trees (Lower ILS): `10-taxon/lower-ILS/estimated-supergenetrees/supergenes-<T>/R<X>/bin.<Z>.txt/RAxML_bootstrap.all` - where X is the replicate number (1...20), T is binning threshold (50 or 75), and Z is the bin id.


### [15-taxon datasets](../binning-response/)

### [Avian and Mammalian datasets (Siavash et al., Science, 2014):](https://www.ideals.illinois.edu/handle/2142/55319) Avian and mammalian simulated datasets, except for the 250bp mammalian model condition.

### [Mammalian (250bp model condition)](https://figshare.com/articles/dataset/plos_datasets/1411146?file=3331334): 20 replicates, 200 genes, 250bp
File contents:

* Gene alignments: `mammalian-250bp/estimated-genetrees/R<X>/<Y>/<Y>.fasta` - where X is the replicate number (1...20), and Y is the gene id.
* Estimated BestML gene trees: `mammalian-250bp/estimated-genetrees/R<X>/<Y>/RAxML_bipartitions.final.f200` - where X is the replicate number (1...20), and Y is the gene id.
* Estimated bootstrap gene trees: `mammalian-250bp/estimated-genetrees/R<X>/<Y>/RAxML_bootstrap.all` - where X is the replicate number (1...20), and Y is the gene id.
* Estimated BestML supergene trees: `mammalian-250bp/estimated-supergenetrees/R<X>/<Y>/RAxML_bipartitions.final.f200` - where X is the replicate number (1...20), and Y is the gene id.
* Estimated bootstrap supergene trees: `mammalian-250bp/estimated-supergenetrees/R<X>/<Y>/RAxML_bootstrap.all` - where X is the replicate number (1...20), and Y is the gene id.
