"bmvnorm.proposal" <-
function(m,n,Sigma) {

  set <- function(Sigma) {
    Sigma <- array(Sigma,c(n,n,m))
    for(k in 1:m) Sigma[,,k] <- chol(Sigma[,,k])
    L <<- matrix(aperm(Sigma,c(1,3,2)),n,m*n)
  }

  get <- function() {
    Sigma <- aperm(array(L,c(n,m,n)),c(1,3,2))
    for(k in 1:m) Sigma[,,k] <- crossprod(Sigma[,,k])
    Sigma
  }
  
  tune <- function(x,scale=1,eps=1.0E-6) {
    Sigma <- array(0,c(n,n,m))
    for(k in 1:m) {
      V <- cov(t(x[k,,]))
      diag(V) <- pmax(diag(V),eps)
      Sigma[,,k] <- scale^2*V
    }
    set(Sigma)
  }

  proposal <- function(mu) {
    mu + matrix(colSums(rnorm(m*n)*L),m,n)
  }

  L <- matrix(0,n,m*n)
  set(Sigma)
  
  list(m=m,
       n=n,
       set=set,
       get=get,
       tune=tune,
       proposal=proposal)
}

