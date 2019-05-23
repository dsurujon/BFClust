addpath('BF_clustering');

% input sequences for existing clustering
treeseqsfile = 'dataset-010-0.fasta';
% new input sequences - to be assigned clusters
newseqsfile = 'dataset-010-1.fasta';
% existing cluster file
clusterdatafile = 'testout/dataset-010-0.mat';
% using all 7 downstream methods 
allmethods = true;

add_to_clustering(newseqsfile, treeseqsfile, clusterdatafile, allmethods);
