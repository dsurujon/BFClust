addpath('BF_clustering');

% input sequences
myfastafile = 'dataset-010-0.fasta';
% number of clusters 
% default is 3000:200:10000(range from 3000 to 10000 for pan-genomes
% with increments of 200). 
% The sequences will be clustered into 3000, 3200, 
% 3400... etc clusters, and the number of clusters that works best 
% will be selected.

% for the sample set, there are 10 clusters so we will scan 5:15 
scan_clusters = 5:1:15;
% number of trees in the forest (10 is recommended)
n_trees = 10;
% output directory
myoutputdir = 'testout/';

% if output directory doesn't exist, create it
if ~exist(myoutputdir, 'dir')
	mkdir(myoutputdir)
end	

% whether to parallelize the BF construction and DM generation steps
% it is HIGHLY recommended to set this to true for improving runtime. 
% When true, n_trees cores will be utilized. 
parallelize = false;

% run the BF clustering pipeline
run_BF_all(myfastafile, scan_clusters , n_trees, myoutputdir, parallelize);