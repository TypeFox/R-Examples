calcFx <- function(pval.mat, x = c(seq(0.001,0.009,by=0.001),seq(0.010,0.985,by=0.005),seq(0.990,0.999,by=0.001))) {
# this function allows us to compute Fxi from the matrix of p-values obtained by many.pval()

  ## retrieve informations from matrix.pval$pvals computed by many.pval()
  nbstats <- ncol(pval.mat)
  N <- nrow(pval.mat)
  law <- unique(rownames(pval.mat))
  if (length(law) > 1) stop("The row names of 'pval.mat' should all be the same!")
  statnames <- colnames(pval.mat)
  
  Fx.mat <- matrix(0,nrow=length(x),ncol=nbstats) 
  for (i in 1:nbstats) {
    Fx.mat[,i] <- (.C("calfx",as.double(pval.mat[,i]),as.integer(N),as.double(x),
                     as.integer(length(x)),fx=as.double(rep(0,length(x)))))$fx
  }
  
  return(structure(list(Fx.mat=Fx.mat,x=x,law=law,statnames=statnames,N=N), class = c("Fx","pvalue","discrepancy","sizepower")))
}

