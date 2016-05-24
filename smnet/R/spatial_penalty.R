
spatial_penalty<-function(adjacency, wgts, lambda, n.segments){
  adj.spam <- make_spam(adjacency)
  pseudo.inds  <- which(colSums.spam(adj.spam) == 1)
  ij.nzero.adj <- triplet(adj.spam)$indices
  in.pseudo    <- ij.nzero.adj[,2] %in% pseudo.inds
  ij.confl     <- ij.nzero.adj[!in.pseudo,]
  n.nzero      <- nrow(ij.confl)
  p.row.ind    <- rep(1:n.nzero, each = 2)
  p.col.ind    <- c(t(ij.confl))
  p.val        <- wgts[rep(ij.confl[,1], each = 2)]*rep(c(-1, 1), n.nzero)
  D2           <- spam(list(i=p.row.ind, j=p.col.ind, p.val), nrow = n.nzero, ncol = n.segments)
  D2           <- t(D2)%*%D2 
  
  if(!is.null(pseudo.inds)){
    ij.pseudo    <- ij.nzero.adj[in.pseudo,]
    n.nzero      <- nrow(ij.pseudo)
    if(is.null(n.nzero)) n.nzero <- 1
    p.row.ind    <- rep(1:n.nzero, each = 2)
    p.col.ind    <- c(t(ij.pseudo))
    if(is.matrix(ij.pseudo)){
      p.val        <- wgts[rep(ij.pseudo[,1], each = 2)]*rep(c(-1, 1), n.nzero)
    }
    if(is.vector(ij.pseudo)){
      p.val        <- wgts[rep(ij.pseudo[1], each = 2)]*rep(c(-1, 1), n.nzero)
    }
    D1           <- spam(list(i=p.row.ind, j=p.col.ind, p.val), nrow = n.nzero, ncol = n.segments)
    D1            <- t(D1)%*%D1
  }  
  return((lambda)*(D2 + D1))  
}
