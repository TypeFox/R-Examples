# See Schilling & Bock (2005, pp. 8-9)

library(testthat)
library(rpf)

context("rescale")

genOrthogonal<-function(dim) { 
  Q<-MOrthogonal(runif(dim))
  return(Q)
}

# Construct an orthogonal matrix whose first few columns are standardized 'M'
# where columns of 'M' are orthogonal.
# Here "standardized 'M'" means each its columns has length 1.
MOrthogonal<-function(M)
{
  # can set the parameter "tol" of "qr" to decide how small value should be 0
  tmp<-qr(M)
  Q<-qr.Q(tmp,complete=TRUE)
  if(is.vector(M)) { if(Q[1]*M[1]<0) Q<- -Q }
  else { if(Q[1,1]*M[1,1]<0) Q<- - Q }
  return(Q)
}

# adapted from clusterGeneration 1.3.1 by Weiliang Qiu, Harry Joe
genPositiveDefMat <- function(dim, low=-1.4, upp=1.4) {
  u<-matrix(0, dim,dim)
  egvalues <- exp(runif(dim,min=low,max=upp))
  diag(u)<-egvalues #the diagonal elements of u are positive
  Sigma<-u
  if(dim>1)
  { Q<-genOrthogonal(dim) # generate an orthogonal matrix 
    Sigma<-Q%*%u%*%t(Q) # the final positive definite matrix
  }
  Sigma
}

for (dims in 1:3) {
  test_that(paste(dims, "dims"), {
    spec <- list()
    spec[[1]] <- rpf.drm(factors=dims, multidimensional=TRUE)
    spec[[2]] <- rpf.grm(factors=dims, outcomes = 3, multidimensional=TRUE)
    spec[[3]] <- rpf.nrm(factors=dims, T.a="random", T.c="random")
    numItems <- length(spec)
    
    test.point <- rnorm(dims)
    param <- list()
    prob <- list()
    for (ix in 1:numItems) {
      param[[ix]] <- rpf.rparam(spec[[ix]])
      prob[[ix]] <- rpf.prob(spec[[ix]], param[[ix]], test.point)
    }
    
    for (ix in 1:numItems) {
      info <- paste("  (While testing item model", class(spec[[ix]]), ")")
      
      nomove <- rep(0, dims)
      padj <- rpf.rescale(spec[[ix]], param[[ix]], nomove, diag(dims))
      prob.adj <- rpf.prob(spec[[ix]], padj, test.point)
      expect_equal(prob.adj, prob[[ix]], 1e-3, label="Unmoved params",
                   info=info)
      
      move <- rnorm(dims)
      padj <- rpf.rescale(spec[[ix]], param[[ix]], move, diag(dims))
      prob.adj <- rpf.prob(spec[[ix]], padj, test.point-move)
      expect_equal(prob.adj, prob[[ix]], 1e-3, label="Moved params", info=info)
      
      cov <- genPositiveDefMat(dims) * lower.tri(diag(dims), TRUE)
      Icov <- t(solve(cov))
      
      padj <- rpf.rescale(spec[[ix]], param[[ix]], nomove, cov)
      prob.adj <- rpf.prob(spec[[ix]], padj, t(test.point %*% Icov))
      expect_equal(prob.adj, prob[[ix]], 1e-3, label="Covariance scaled params", info=info)
      
      padj <- rpf.rescale(spec[[ix]], param[[ix]], move, cov)
      prob.adj <- rpf.prob(spec[[ix]], padj, t(t(test.point-move) %*% Icov))
      expect_equal(prob.adj, prob[[ix]], 1e-3,
                   label="Moved and covariance scaled params", info=info)
    }
  })
}

if (0) {
  sum((param[[ix]][1:2]) * test.point)
  sum((param[[ix]][1:2] %*% cov) * (test.point %*% t(solve(cov))))
}
