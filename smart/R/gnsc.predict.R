gnsc.predict <-
function(x, x.class=NULL, all.mean, tilde.mu, row.struc, col.struc, inv.cov.mat, path.se.id){
	x=as.matrix(x)
	d=ncol(x)
	class = unique(col.struc)
	M = length(class)
	path = unique(row.struc)
	K = length(path)
	n = length(col.struc)
  clust.n = gnsc.restruc(col.struc)$path.n

	e.class = rep(0, d)
	x2=x-all.mean
	score=matrix(0,nrow=d,ncol=M)
	for(ci in 1:M){
		x3=x2-tilde.mu[,ci]
		G.norms = matrix(0,nrow=K,ncol=d)
		for(k in 1:K){
			tmp.id = (path.se.id[k,1]:path.se.id[k,2])
			G.x3 = x3[tmp.id,]
			if(!is.matrix(G.x3)) G.x3 = matrix(G.x3,1,length(G.x3)) 
			G.norms[k,] = diag(t(G.x3)%*%inv.cov.mat[[k]]%*%G.x3)
		}
		score[,ci]=colSums(G.norms)-2*log(clust.n[ci]/n)
	}
	e.class = apply(score, 1, minid)
	
	if(is.null(x.class)) error=NULL
	if(!is.null(x.class)) error=sum(class[e.class]!=x.class)/d
	output = list()
	output$yhat = e.class
	output$error = error
	return(output)
}
