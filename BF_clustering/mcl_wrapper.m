function [gnew, msg] = mcl_wrapper(dm)
	W = exp(-dm.^2)>0.9;
	adj = W;

    [gnew, msg] = mcl_ds(adj,0,0,0,0,50);

end