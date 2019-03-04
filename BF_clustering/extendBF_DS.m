function all_clusters = extendBF_DS(clusters_tree,tree,treerefs)
    all_clusters = zeros(1,length(treerefs));
    for i=1:length(clusters_tree)
        tree_ix=tree{i}{1};
        all_clusters(treerefs==tree_ix)=clusters_tree(i);
    end

end