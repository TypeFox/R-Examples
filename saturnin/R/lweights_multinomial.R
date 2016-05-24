lweights_multinomial <-
function(data,prior = defaut.prior,nbcores = 1){
  p <- ncol(data)
  if (missing(prior)){
    r <- max(data)
    defaut.prior <- prior_unif_dirichlet(p,r,0.5*r^2 )
  } else {
    r <- dim(prior)[1] 
  }
  
  dataf <- as.data.frame(data)
  for (i in 1:p){
    dataf[,i] <- factor(dataf[,i],levels=1:r,ordered = TRUE)
  }
  
  uptri <- upper.tri(matrix(0,p,p))
  
  if (requireNamespace("parallel",quietly = TRUE) &&
        (nbcores > 1)){
    lW <- matrix(parallel::mcmapply(function(y, i, j) if (y){
      post.ij <- table(dataf[,i],dataf[,j]) + prior[,,i,j]
      sum.intdiff.lgamma(prior[,,i,j],post.ij)    
    } else 0,
    uptri, 
    row(uptri), 
    col(uptri),
    mc.cores = nbcores), 
    nrow = nrow(uptri))  
  } else {
    lW <- matrix(mapply(function(y, i, j) if (y){
      post.ij <- table(dataf[,i],dataf[,j]) + prior[,,i,j]
      sum.intdiff.lgamma(prior[,,i,j],post.ij)    
    } else 0,
    uptri, 
    row(uptri), 
    col(uptri)), 
    nrow = nrow(uptri))
  }
  
  diaglW <- sapply(1:p,function(i) sum.intdiff.lgamma(diag(prior[,,i,i]),diag(prior[,,i,i])+table(dataf[,i])))
  lW <- lW + t(lW)
  lW <- sapply(1:p,function(x) lW[,x] - diaglW[x])
  lW <- t(sapply(1:p,function(x) lW[x,] - diaglW[x]))
  diag(lW) <- 0
  return(lW)
}
