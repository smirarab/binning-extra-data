# Simulation Details

## 0- Species tree

The first task that needs to be done is generating a species tree with branch lengths in coalescence units (for internal branches). To do this, we use our biological datasets and draw coalescent unit branch lengths on a fixed topology using estimated gene trees. We use MP-EST to draw branch length.

For the mammalian simulation, we used the biological MP-EST tree (unbinned) as the model tree. For simulating the avian dataset, we used the binned MP-EST tree on the TENT matrix from the avian project but re-estimated the branch lengths on that model tree using only the longest genes with at least 10,000 sites. This was necessary because most gene trees had very low support values (and thus high estimation error); branch lengths estimated in coalescent units are directly impacted by observed gene tree discordance, and including gene trees with high estimation error (i.e. low support gene trees) inflates the amount of ILS.

The model species trees for 1X model condition: [avian](avian-model-species.tre) and [mammalian](mammalian-model-species.tre).

## 1- Simulating Gene trees

We ran

    python generateCoalescentTrees.py inputSpeciesTreeWithBranchLengthInCoalescentUnits 1000 simulatedGeneTrees 
to simulate 1000 gene trees down the species given by `inputSpeciesTreeWithBranchLengthInCoalescentUnits`. The script [generateCoalescentTrees.py](generateCoalescentTrees.py) uses the command `treesim.contained_coalescent` from the dendropy package to simulate gene trees according to the species tree with coalescent unit branch lengths. Running this script requires that DendroPy is available on the system.

## 2- Adjusting branch lengths

Now we change the branch lengths into expected numbers of substitutions. We run:

    ./transposeBranchLengthsByPercentile input.trees.file=biologicalGeneTrees  input.trees.to.change.file=simulatedGeneTrees  output.file=simulatedTrees_realisticBranchLengths 
where biologicalGeneTrees contains the real gene trees estimated from actual alignments. For the avian data set, we used the gene trees reconstructed from the 190 longest introns. For the mammalian dataset, we used all 447 gene trees.

[transposeBranchLengthsByPercentile](transposeBranchLengthsByPercentile) is a program that requires the [Bio++ package](http://biopp.univ-montp2.fr/) to compile (at the beginning of the C++ code, there is a command line for compiling it).

## 3- Simulating gene sequences

To simulate gene sequences, we used [bppseqgen](http://home.gna.org/bppsuite/doc/bppseqgen.html), a member of [bppsuite](http://home.gna.org/bppsuite/).

To simulate dna sequences according to parameters estimated from the avian data, we used:

    bppseqgen number_of_sites=1500 input.tree.file=geneTree param=bppml.estimates output.sequence.file=sequences.fasta 
where bppml.estimates is a file containing the following parameters (estimated from the avian data and used for both avian and mammalian datasets):

	# Substitution model parameters:
	model = GTR(a=1.062409952497, b=0.133307705766, c=0.195517800882, d=0.223514845018, e=0.294405416545, theta=0.469075709819, theta1=0.558949940165, theta2=0.488093447144)
	
	# Rate distribution parameters:
	rate_distribution = Gamma(n=4, alpha=0.370209777709)
This file is returned by the program `bppml`, also from `bppsuite`. We ran bbpml on the subset of 1185 avian genes that included all the taxa.