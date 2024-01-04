Dataset
===

The following files are provided:

* [MCcoal.zip](MCcoal.zip) includes the MCcoal control files and the seeds used 

* Alignments, true gene trees, and estimated gene trees for 10-taxon and 15-taxon datasets (5-taxon was simulated by L&E and provided to us): 
	* Gene alignments: [t10_t15_gene_alignments.tar.bz2](t10_t15_gene_alignments.tar.bz2)
	* True gene trees: [t10_t15_true_gene_trees.tar.bz2](t10_t15_true_gene_trees.tar.bz2) 
	* Estimated gene trees: [t10_t15_estimated_genetrees.tar.bz2](t10_t15_estimated_genetrees.tar.bz2)

* Supergene alignments and trees:
	* Supergene alignments: for each supergene, the ailgnment and the partition file is give: [supergene_alignments.tar.bz2](supergene_alignments.tar.bz2)
	* Supergene trees: [supergene_trees.tar.bz2](supergene_trees.tar.bz2)

* MP-EST results (including the input used to each MP-EST run): [MPEST_input_and_output.tar.bz2](MPEST_input_and_output.tar.bz2)


Simulation Procedure
===

Gene trees were simulated using MCcoal, with control files given in our dataset. These control files define the species tree. The species trees are in the caterpillar form (shown below):

* 10 species:

```
(((((((((A #.05,B #.05):0.005 #.05,C #.05):0.01 #.05, D #.05):0.015 #.05,E #.05):0.02 #.05,F #.05):0.025 #.05,G #.05):0.03 #.05,H #.05):0.035 #.05,I #.05):0.04 #.05, J #.05):0.54 #.05;
```

* 15 species:

```
((((((((((((((A #.05,B #.05):0.005 #.05,C #.05):0.01 #.05, D #.05):0.015 #.05,E #.05):0.02 #.05,F #.05):0.025 #.05,G #.05):0.03 #.05,H #.05):0.035 #.05,I #.05):0.04 #.05,J #.05):0.045 #.05,K #.05):0.05 #.05,L #.05):0.055 #.05,M #.05):0.06 #.05,N #.05):0.065 #.05,O #.05):0.565 #.05;
```

To run MCcoal, the following command was used:

```
   printf "10000 1000" | PATH_TO_MCCOAL/MCcoal
```

This simulated 10,000 gene trees, which we divided into 10 replicates of 1000 genes or 100 genes. 

For each true gene tree, we then simulated alignments using `bppseqgen`, using the following command:

```    
mkdir allTrees
split -a 4 -l 1 out.trees 
for i in x* ; do mv $i  allTrees/ ; done
for i in allTrees/x* ; do bppseqgen number_of_sites=1000 input.tree.file=$i param=bpp.options output.sequence.file=$i".fasta" ; done
```

The file `bpp.options` is the same as what was used in our statistical binning paper: 

```   
# Substitution model parameters:
model = GTR(a=1.062409952497, b=0.133307705766, c=0.195517800882, d=0.223514845018, e=0.294405416545, theta=0.469075709819, theta1=0.558949940165, theta2=0.488093447144)

# Rate distribution parameters:
rate_distribution = Gamma(n=4, alpha=0.370209777709)
```


Estimating gene trees
===
To estimate gene trees, we used RAxML, version 8.0.19. We used the following commands. 

* Unbinned gene trees (unpartitioned analyses):
	* maximum likelihood analyses:
	
          raxmlHPC-8.0.19-SSE3 -m GTRGAMMA -n best -s [alignment_file] -N 10 -p [random_seed_number]
    
	* bootstrapping

          raxmlHPC-8.0.19-SSE3 -m GTRGAMMA -n ml -s [alignment_file] -N 100 -b [random_seed_number] -p [random_seed_number]

	* drawing support values onto the maximum likelihood tree
	
          raxmlHPC-8.0.19-SSE3 -f b -m GTRGAMMA -n final.f100 -z RAxML_bootstrap.ml -t RAxML_bestTree.best


* Supergene trees (partitioned analyses):

	* bootstrapping
	
          raxmlHPC-8.0.19-SSE3 -m GTRGAMMA -n ml -s [alignment_file] -N 100 -b [random_seed_number] -p [random_seed_number] -M -q supergene.part


Note that the partition files are provided as part of our dataset. 


Binning
===
To perform binning, we used a pipeline available on [gitub](https://github.com/smirarab/binning). As input to the binning pipeline, we used the `RAxML_bipartitions.final` files produced as part of our estimation of unbinned gene trees. We varied the threshold (25, 50, and 75). 

For the 75% threshold, we also ran some tests where we investigated the impact of breaking ties in the binning code in multiple ways. 
To force the binning pipeline to break ties differently from one run to another, we modified the pipeline to create random orderings of gene names that are used as input to the pipeline.
Thus, instead of 

    ls| grep -v ge|sed -e "s/.50$//g"> genes

we used:

    ls| grep -v ge|sed -e "s/.50$//g"|sort -R > genes


Species tree estimation
===
To estimate species trees, we followed the following steps. 

1. We rooted gene trees (using custom scripts based on dendropy) on the outgroup (E for 5-taxon, J for 10-taxon, and O for 15-taxon)
2. We created 100 bootstrap replicate inputs for MP-EST. To do this, we matched the first bootstrap replicate of each input (super-) gene tree to build input 1 (BS.1), matched second bootstrap replicate to create input 2 (BS.2), and so on. These files (BS.1, BS.2, ..., BS.100) are provided as part of our dataset. 
3. We ran MP-EST on each input (BS.i). We ran MP-EST 10 times on each input an picked the result with the highest likelihood value. The final tree from the best run is provided in our datasets
4. We estimated the greedy consensus of all 100 bootstrap runs, and used this consensus as the estimate of the species tree. 

