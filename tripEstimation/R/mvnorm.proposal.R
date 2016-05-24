"mvnorm.proposal" <-
function(m,n,Sigma) {
  set <- function(Sigma) {
    L <<- chol(Sigma)
  }

  get <- function() {
    crossprod(L)
  }
  
  tune <- function(x,scale=1,eps=1.0E-6) {
    d <- dim(x)
    Sigma <- cov(t(matrix(x,d[1]*d[2],d[3])))
    diag(Sigma) <- pmax(diag(Sigma),eps)
    set(scale^2*Sigma)
  }

  proposal <- function(mu) {
    mu+matrix(rnorm(n)%*%L,m,n)
  }

  L <- matrix(0,m*n,m*n)
  set(Sigma)
  
  list(m=m,
       n=n,
       set=set,
       get=get,
       tune=tune,
       proposal=proposal)
}

