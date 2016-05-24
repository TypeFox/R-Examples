"norm.proposal" <-
function(m,n,sigma) {

  set <- function(sigma) {
    s <<- sigma
  }
  
  get <- function() {
    s
  }

  set <- function(sigma) {
    s <<- sigma
  }
  
  tune <- function(x,scale=1,eps=1.0E-6) {
    set(scale*pmax(apply(x,1:2,sd),eps))
  }

  proposal <- function(mu) {
    matrix(rnorm(m*n,mu,s),m,n)
  }

  s <- sigma
  
  list(m=m,
       n=n,
       set=set,
       get=get,
       tune=tune,
       proposal=proposal)
}

