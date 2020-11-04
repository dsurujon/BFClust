# BFClust

Boundary-Forest Clustering is a pan-genome clustering pipeline written in MATLAB. Boundary-Forest Clustering is done in 3 major steps:    
1. **Boundary-Forest.** This generates a number of Boundary-Trees based on sequence similarity (See [Mathy et al., 2015](https://arxiv.org/abs/1505.02867) for more detail on Boundary-Forest). A sequence is either added to the tree as a new representative, or if there is an existing representative on the tree that is sufficiently close to this sequence, the sequence is annotated with this representative, and omitted from the tree. Boundary-Trees thus contain a small subset of input sequences as representative sequences, that are arranged in a tree-structure based on sequence similarity.
2. **Cluster.** Clustering is performed on each Boundary-Tree. Currently 7 clustering methods are implemented downstream of Boundary-Forest. These include hierarchical, kmeans, kmeans(vectorized), spectral, spectral (Shi-Malik normalized), spectral (Ng-Jordan-Weiss normalized) and Markovian clustering (MCL). We recommend using MCL as the main method. Once the representative sequences are clustered, the clustering assignments are extended to the full dataset. 
3. **Consensus Clustering.** The clustering on the different Boundary-Tree representatives may not always yield identical results. Taking a consensus of the clustering assignments across the forest reduces errors for all downstream clustering methods. A consensus score for each element and each cluster will also be generated as a quality metric. The consensus score for a cluster or an item is a value between 0 and 1. A score of 1 indicates perfect agreement among individual clustering assignments.      

The major advantage of BFClust is that it outputs the level of certainty (a consensus score) associated with each item and/or each cluster. This gives a measure of cluster "quality" when we do not know what the "real" clusters are supposed to be. Another advantage of BFClust is that it stores the Boundary-Forest, making it possible to add new sequences to the clustering without having to alter the existing clustering assignments ("cluster augmentation"). This not only reduces the time necessary to obtain cluster assignments for an incoming set of sequences (e.g. a newly sequenced bacterial isolate), but also keeps the existing clustering assignments the same.    
    
There are three main scripts that can be used. ```run_BF_all.m``` and ```run_BF_single.m``` are for clustering a new dataset *de novo*, and ```add_to_clustering.m``` is for adding new sequences to an existing clustering partition. 

## *De novo* clustering
This is used for clustering a new sequence set. BFClust has the option of clustering a sequence set using one (or all) of 7 clustering algorithms: hierarchical, kmeans, kmeans(vectorized), spectral, spectral (Shi-Malik normalized), spectral (Ng-Jordan-Weiss normalized) and Markovian clustering. ```run_BF_all.m``` runs the entire clustering pipeline, and all 7 clustering algorithms. On the other hand, ```run_BF_single.m``` runs only one specified clustering algorithm. The input arguments for these two scripts are as follows:    
* ```fastafile```: name of the fasta file containing all amino acid sequences to be clustered
* ```krange```: list of values for k, the number of clusters, to be scanned
* ```replicates```: size of the boundary forest (i.e. how many trees will be generated)
* ```outdir```: name of the output directory
* ```isparallel```: whether ot not to parallelize the boundary forest construction, and distance matrix calculation. It is HIGHLY recommended that this is set to ```true``` to reduce runtime. When set to ```true```, the number of cores requires will be equal to ```replicates``` 
* ```methodname```: **(Only for ```run_BF_single```)** Name of clustering method to be used (one of ```HIE, KMN, KMV, SP1, SP2, SP3, MCL```.     
    
    
There will be up to 10 output files generated inside ```outdir```, named after the input file name. The ```[datasetname].mat``` file contains all intermediate data, including boundary forest, distance matrices, clustering assignments, consensus.     
The ```[datasetname].csv``` file is a table of sequence headers (from the input fasta file), the cluster assignments, and the consensus score for each item and its cluster assignment. When ```cluster_BF_all``` is used, there will be 7 columns for cluster assignments and 7 columns for consensus score, each corresponding to one method.     
The ```[datasetname][Methodname].csv``` files contain the cluster consensus scores for each method.     
Finally, the ```.png``` file shows the SSE trace for each tree for the selected method(s) - with the exception of MCL. For MCL, the best number of clusters is determined without the use of elbow detection on the SSE traces, therefore no plot is generated.     
An example use case is provided in the script ```cluster_example.m```. This example takes the  ```dataset-010-0.fasta``` sequence set as input (there are 500 sequences; 50 copies of 10 genes with a small mutation rate, so this should yield 10 clusters). The example output of BFClust can be found in ```testout```. 


## Cluster augmentation (Adding to existing clustering)
This is used when a clustering partition already exists, and one wishes to assign clusters to a new sequence set. This is especially useful when a large number of sequences have already been clustered, and a relatively small sequence set is to be assigned clusters. The advantage here is three-fold:     
1. Existing cluster assignments are not changed
2. Adding new sequences is faster than clustering the old and new sequences together    
3. The confidence scores are computed again, giving a level of uncertainty for the newly added clustering      
The ```add_to_clustering.m``` script assigns clusters to a new set of sequences by finding the closest cluster in the existing clustering of the existing set of sequences.     
The inputs this function takes are as follows:    
* ```newseqsfile```: name of fasta file containing new sequences
* ```treeseqsfile```: name of fasta file containing the sequences that have already been clustered
* ```clusterdatafile```: name of ```.mat``` file containing the clustering data for the existing set of sequences (this file is generated during *de novo* clustering)
* ```allmethods```: whether or not all 7 clustering methods are considered (i.e. is ```clusterdatafile``` generated using ```cluster_BF_all``` or ```cluster_BF_single```?)    
    
An example use case is provided in the script ```add_to_cluster_example.m```. This example takes the  ```dataset-010-1.fasta``` sequence set as new input (similar to ```dataset-010-0.fasta```, there are 500 sequences; 50 copies of 10 genes with a small mutation rate), and adds it to the clustering results of ```dataset-010-0```. The example output is also in ```testout```. 
    
## Citing
Please cite the BFClust preprint:     
    
[Boundary-Forest Clustering: Large-Scale Consensus Clustering of Biological Sequences](https://www.biorxiv.org/content/10.1101/2020.04.28.065870v1)
Defne Surujon, José Bento, Tim van Opijnen
bioRxiv 2020.04.28.065870; doi: https://doi.org/10.1101/2020.04.28.065870

-----------------
Several scripts utilize existing code from others' libraries    
* Spectral Clustering is modified from https://www.mathworks.com/matlabcentral/fileexchange/34412-fast-and-efficient-spectral-clustering
* MCL is modified from https://github.com/biocoder/SBEToolbox/blob/master/mcl.m
* cell2csv: https://www.mathworks.com/matlabcentral/fileexchange/7601-cell2csv
* Hungarian matching algorithm (lapjv.m): https://www.mathworks.com/matlabcentral/fileexchange/26836-lapjv-jonker-volgenant-algorithm-for-linear-assignment-problem-v3-0
