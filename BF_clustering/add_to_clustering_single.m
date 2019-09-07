% DS Feb 2019
% This function takes in a new fasta file, and an existing clustering - single downstream method
% partition (a .mat data file, output of runBF_all), and assigns the best
% cluster to each of the incoming sequences. 
% [input] newseqsfile: fasta file of new sequences
% [input] treeseqsfile: fasta file used for the existing clustering
% [input] clusterdatafile: .mat file containing clustering information for
% sequenes in treeseqsfile. This should at least have the variables
% 'trees', 'clusterres_ext' and 'consclust'
% [output] new_consclust: final consensus clustering of new sequences. 
% [output] new_clusterres_ext: best cluster assignment for each tree (this
% is used for consensus clustering to get to new_consclust). 

function [new_consclust, new_clusterres_ext] = add_to_clustering_single(newseqsfile, treeseqsfile, clusterdatafile)

newseqs = read_and_clean(newseqsfile);
% make a boundary tree out of newseqsfile
eps = 0.1;
max_deg = 10;
[newtree,newtree_node_ref,newdata_order_ix] = boundary_tree(newseqs, eps, max_deg);
% grab the sequences present in newtree
newtreeseqs = newseqs(arrayfun(@(x) newtree{x}{1}, [1:length(newtree)]));


treeseqs = read_and_clean(treeseqsfile);
load(clusterdatafile, 'trees', 'clusterres_ext', 'consclust');
ntree = size(trees, 2);

% assign to a representative for each BT

new_clusterres_ext = zeros(length(newtreeseqs), ntree);
for tree_ix = 1:ntree 
	[rep_seq, rep_seq_tree_ix] = find_closest_on_BF(newtreeseqs,trees{tree_ix},treeseqs);
	clusters_thistree = clusterres_ext(rep_seq, tree_ix);
	new_clusterres_ext(:, tree_ix)=clusters_thistree;
	
end
% find closest neighbor on consensus. 
Idx = knnsearch(clusterres_ext, new_clusterres_ext);
% grab the consensus cluster assignment
cons_assn = consclust(Idx);
%extend consensus clustering
cons_assn_ext = extendBF_DS(cons_assn, newtree, newtree_node_ref);
new_consclust = cons_assn_ext;





end