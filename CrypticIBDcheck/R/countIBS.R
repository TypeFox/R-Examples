countIBS <- function(x) {
  if (!is(x, "snp.matrix"))
    stop("argument must be a snp.matrix object")
  n <- nrow(x)
  mni0=matrix(0,n,n) #matrix of count of snps with IBS=0 for pairs of individuals
  mni1=matrix(0,n,n)  
  mni2=matrix(0,n,n)
  out <- .C("countIBS", t(x@.Data), as.integer(n),  as.integer(ncol(x)),
          as.integer(mni0), as.integer(mni1), as.integer(mni2), PACKAGE="CrypticIBDcheck")
  out[[4]] <- matrix(out[[4]], nrow=n, byrow=TRUE)
  out[[5]] <- matrix(out[[5]], nrow=n, byrow=TRUE)
  out[[6]] <- matrix(out[[6]], nrow=n, byrow=TRUE)
  names(out)<-c(rep(NA,3),"IBS0","IBS1","IBS2")

  out[4:6]
}

