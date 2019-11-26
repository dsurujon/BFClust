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

function[item_consensus, cluster_consensus] = consensus_quality_single(clusterres_ext, consclust, n_tree)
    N = size(clusterres_ext,1);
    
    %make the connectivity matrices
    M = zeros(N, N, n_tree);
    for i=1:n_tree
        M(:,:,i) = clusterres_ext(:,i)==clusterres_ext(:,i)';
    end
    % make consensus - average of all the M
    M_cons = mean(M,3);
    
    % item consensus
    item_consensus = zeros(1,N);
    for i = 1:N
       thiscluster = consclust(i);
       ix = consclust==thiscluster;
       M_cons_sub = M_cons(i, ix);
       item_consensus(i) = mean(M_cons_sub);
    end
    
    % cluster consensus
    allclusters = unique(consclust)';
    cluster_consensus = zeros(1,length(allclusters));
    i=1;
    for cluster=allclusters
        ix = consclust==cluster;
        M_cons_sub = M_cons(ix,ix);
        cluster_consensus(i) = mean(mean(M_cons_sub));
        i=i+1;
    end

end