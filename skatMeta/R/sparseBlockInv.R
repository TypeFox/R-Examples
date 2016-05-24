# sparseBlockInv <- function(Theta,thresh = 2*2^{-6}){
	# require(Matrix)
	# require(igraph)
	# G <- graph.adjacency(as.matrix(Theta) > thresh,mode = "undirected")
	# g <- igraph:::clusters(G)
	# Theta.i.list <- as.list(by(1:nrow(Theta),g$membership,function(inds){
		# as.matrix(solve(Theta[inds,inds]))
	# }) )
	# Theta.i.perm <- bdiag(Theta.i.list)
	# perm <- match(1:nrow(Theta),unlist(tapply(1:nrow(Theta),g$membership,function(x){x})))
	# return(Theta.i.perm[perm,perm])
# }