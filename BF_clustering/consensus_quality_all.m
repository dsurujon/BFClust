% Consensus quality 
% Compute a consensus score for each item and/or each cluster, indicating
% how much agreement there is among the clustering partitions we get from
% individual trees. 
% For each partition, we compute a connectivity matrix M, Where M(i,j)=1 if
% i and j belong to the same cluster, 0 otherwise. Then we take the average
% of the connectivity matrices to obtain M_cons, the consensus matrix. 
% Cluster consensus is the average of M_cons(ix,ix) where ix are the
% indices of items belonging to this cluster. 
% Item consensus is the average of M_cons(i, ix) where i is the index of
% the item, and ix are the indices of all other items that belong to the
% same consensus cluster as item i. 
% [input] clusterres_ext: 1x7 cell where each cell has the clustering
% output of one clustering method (NxN_tree numeric matrix)
% [input] consclust: 1x7 cell that has consensus clustering result for 
% each of the 7 methods
% [output] item_consensus: 1x7 cell (one for each method)
% [output] cluster_consensus: 1x7 cell (one for each method)

function[item_consensus, cluster_consensus] = consensus_quality_all(clusterres_ext, consclust, n_tree)
N = size(clusterres_ext{1},1);

item_consensus = cell(1,7);
cluster_consensus = cell(1,7);
for method = 1:7 
    clusterres_method = clusterres_ext{method};
    consensus_method = consclust{method};
    
    allclusters = unique(consensus_method)';

    item_scores_method = zeros(1,N);
    cluster_scores_method = zeros(1,length(allclusters));
    
    i=1;
    for cluster = allclusters
        %find the indices of this cluster
       thiscluster_ix = consensus_method==cluster;
       n_clust = sum(thiscluster_ix);
       % make the connectivity matrices (only for this cluster to save
       % memory)
       M = zeros(n_clust, n_clust, n_tree); 
       for tree=1:n_tree
           M(:,:,tree) = clusterres_method(thiscluster_ix,tree)==clusterres_method(thiscluster_ix,tree)';
       end
       %make the consensus matrix by averaging the connectivity matrices
       M_cons = mean(M,3);

       cluster_scores_method(i) = mean(mean(M_cons));

       item_scores_method(thiscluster_ix) = mean(M_cons);
       i = i+1;
    end

    item_consensus{method} = item_scores_method;
    cluster_consensus{method} = cluster_scores_method;

end

end