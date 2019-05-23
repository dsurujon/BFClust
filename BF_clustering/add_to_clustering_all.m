% DS Feb 2019
% This function takes in a new fasta file, and an existing clustering - 7 downstream methods
% partition (a .mat data file, output of runBF_all), and assigns the best
% cluster to each of the incoming sequences. 
% [input] newseqsfile: fasta file of new sequences
% [input] treeseqsfile: fasta file used for the existing clustering
% [input] clusterdatafile: .mat file containing clustering information for
% sequenes in treeseqsfile. This should at least have the variables
% 'trees', 'clusterres_ext' and 'consclust'
% [output] new_consclust: final consensus clustering of new sequences. Cell
% i has the output for method i (Ward/Kmeans/...etc)
% [output] new_clusterres_ext: best cluster assignment for each tree (this
% is used for consensus clustering to get to new_consclust). Cell
% i has the output for method i (Ward/Kmeans/...etc)

function [new_consclust, new_clusterres_ext] = add_to_clustering_all(newseqsfile, treeseqsfile, clusterdatafile)

newseqs = fastaread(newseqsfile);
% make a boundary tree out of newseqsfile
eps = 0.1;
max_deg = 10;
[newtree,newtree_node_ref,newdata_order_ix] = boundary_tree(newseqs, eps, max_deg);
% grab the sequences present in newtree
newtreeseqs = newseqs(arrayfun(@(x) newtree{x}{1}, [1:length(newtree)]));


treeseqs = fastaread(treeseqsfile);
load(clusterdatafile, 'trees', 'clusterres_ext', 'consclust');
ntree = size(trees, 2);

new_consclust = cell(1,7);
% assign to a representative for each method, each BT
new_clusterres_ext = cell(1,7); 
for method_ix = 1:7
    ext_thismethod = zeros(length(newtreeseqs), ntree);
   for tree_ix = 1:ntree 
        [rep_seq, rep_seq_tree_ix] = find_closest_on_BF(newtreeseqs,trees{tree_ix},treeseqs);
        clusters_thistree = clusterres_ext{method_ix}(rep_seq, tree_ix);
        ext_thismethod(:, tree_ix)=clusters_thistree;
        
   end
   new_clusterres_ext{method_ix} = ext_thismethod;
   % find closest neighbor on consensus. 
    Idx = knnsearch(clusterres_ext{method_ix}, new_clusterres_ext{method_ix});
    % grab the consensus cluster assignment
    cons_assn = consclust{method_ix}(Idx);
    %extend consensus clustering
    cons_assn_ext = extendBF_DS(cons_assn, newtree, newtree_node_ref);
    new_consclust{method_ix} = cons_assn_ext;
end




end