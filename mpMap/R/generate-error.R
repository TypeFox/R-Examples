generate_error <- function(geno, error.prob)
{
  obsgeno <- geno
  n.founders <- nrow(geno$founders)
  n.mrk <- ncol(geno$founders)
  n.finals <- nrow(geno$finals) 

  fdr.err <- matrix(data=sample(c(TRUE,FALSE), n.founders*n.mrk, replace=TRUE, 
	prob=c(error.prob, 1-error.prob)), nrow=n.founders, ncol=n.mrk)
  fin.err <- matrix(data=sample(c(TRUE,FALSE), n.finals*n.mrk, replace=TRUE,
	prob=c(error.prob, 1-error.prob)), nrow=n.finals, ncol=n.mrk)

  # replacing IBD genotypes with something different
  for (i in 1:n.founders) {
     prob <- rep(error.prob/(n.founders-1), n.founders) 
     prob[i] <- 0

     obsgeno$founders[fdr.err==1 & obsgeno$founders==i] <- sample(1:n.founders, sum(fdr.err==1 & obsgeno$founders==i), replace=TRUE, prob=prob)
     obsgeno$finals[fin.err==1 & obsgeno$finals==i] <- sample(1:n.founders, sum(fin.err==1 & obsgeno$finals==i), replace=TRUE, prob=prob)
  }

  return(obsgeno)
}

