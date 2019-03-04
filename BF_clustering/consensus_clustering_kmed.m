% Consensus clustering
%
% given a set of clusterings (as a matrix with each object as a row, and
% each clustering as a column), find a consensus clustering. The cluster
% assignments in clusterres_ext are used as feature vectors and
% clustered using kmedioids into a new set of clusters. 
%
% [input] clusterres_ext: n X m matrix with n data points, and m clustering
% assignments
% [input] bestcluster: 1xm list of number of clusters used for each column
% (clustering) in clusterres_ext
% [output] consclust: consensus clustering assignments

function [consclust, kcons] = consensus_clustering_kmed(clusterres_ext, bestcluster)
kcons = mode(bestcluster);
disp(kcons);
consclust=kmedoids(clusterres_ext,kcons,'Distance','hamming');
end