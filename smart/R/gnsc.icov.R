gnsc.icov <-
function(x,rname=rownames(x), pathway, path.n, path.se.id){
	d = dim(x)
	path = unique(pathway)
	K = length(path)
	cov.mat = list()
	inv.cov.mat=list()
	length(inv.cov.mat)=K
	s0.diag=rep(0,length(pathway))
	for(k in 1:K){
		tmp.id = (path.se.id[k,1]:path.se.id[k,2])		
		tmp.mat = x[tmp.id,]
		if(!is.matrix(tmp.mat)) tmp.mat = matrix(tmp.mat,1,length(tmp.mat))
		cov.mat[[k]]=cov(t(tmp.mat))
		s0.diag[tmp.id] = diag(cov.mat[[k]])	
	}
	s0=median(s0.diag)
	for(k in 1:K){
		inv.cov.mat[[k]]=solve( cov.mat[[k]] + diag(s0, path.n[k]) )
	}
	return(inv.cov.mat)
}
