% DS July 2018
% Run each point of a new sequence set through a boundary tree
% and find its closest representative on the tree

function [X,Y] = extendBF_v3(items, tree, treeseqs)
    n = length(items);
    X = zeros(n,1);
	Y = zeros(n,1);
    % for each new item
    for i=1:n
		%disp(i);
        curr_tree_node_ix = 1;
        curr_node_all_children_node_ix = tree{curr_tree_node_ix}{2};
        
        curr_node_dist = seqpdist([treeseqs(tree{curr_tree_node_ix}{1}),...
            items(i)],'Method','alignment-score',...
            'ScoringMatrix','BLOSUM62','GapOpen',10,'ExtendGap',0.5,...
            'PairwiseAlignment',true);
        
        path_ixs = [curr_tree_node_ix];
        path_dists = [curr_node_dist];
        while ~isempty(curr_node_all_children_node_ix)
            smallest_dist = 1000000;
            for child_node_ix = curr_node_all_children_node_ix  
                dis = seqpdist([treeseqs(tree{child_node_ix}{1}),...
                    items(i)],'Method','alignment-score',...
                    'ScoringMatrix','BLOSUM62','GapOpen',10,'ExtendGap',0.5,...
                    'PairwiseAlignment',true);
                
                if (dis < smallest_dist)
                    smallest_dist = dis;
                    best_child_node_ix = child_node_ix;
                end    
            end
            path_ixs = [path_ixs best_child_node_ix];
            path_dists = [path_dists smallest_dist];
            curr_tree_node_ix = best_child_node_ix;
            curr_node_all_children_node_ix = tree{curr_tree_node_ix}{2};
        end
        
        [minpath, argminpath] = min(path_dists);
        X(i) = tree{path_ixs(argminpath)}{1};
		Y(i) = path_ixs(argminpath);
    end
end