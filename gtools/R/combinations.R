# $Id: combinations.R 1083 2007-03-23 22:53:00Z warnes $
#

##
## From email by Brian D Ripley <ripley@stats.ox.ac.uk> to r-help
## dated Tue, 14 Dec 1999 11:14:04 +0000 (GMT) in response to
## Alex Ahgarin <datamanagement@email.com>.  Original version was
## named "subsets" and was Written by Bill Venables.  
##

combinations <- function(n, r, v = 1:n, set = TRUE, repeats.allowed=FALSE)
{
  if(mode(n) != "numeric" || length(n) != 1 
     || n < 1 || (n %% 1) != 0) stop("bad value of n") 
  if(mode(r) != "numeric" || length(r) != 1 
     || r < 1 || (r %% 1) != 0) stop("bad value of r") 
  if(!is.atomic(v) || length(v) < n) 
    stop("v is either non-atomic or too short")
  if( (r > n) & repeats.allowed==FALSE)
    stop("r > n and repeats.allowed=FALSE")
  if(set) {
    v <- unique(sort(v))
    if (length(v) < n) stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  ## Inner workhorse
  if(repeats.allowed)
    sub <- function(n, r, v)
      { 
        if(r == 0) v0 else
        if(r == 1) matrix(v, n, 1) else
        if(n == 1) matrix(v, 1, r) else
        rbind( cbind(v[1], Recall(n, r-1, v)),
              Recall(n-1, r, v[-1]))
      }
  else
    sub <- function(n, r, v)
      { 
        if(r == 0) v0 else
        if(r == 1) matrix(v, n, 1) else
        if(r == n) matrix(v, 1, n) else
        rbind(cbind(v[1], Recall(n-1, r-1, v[-1])),
              Recall(n-1, r, v[-1]))
      }
  sub(n, r, v[1:n])
}

##
## Original version by Bill Venables and cited by by Matthew
## Wiener (mcw@ln.nimh.nih.gov) in an email to R-help dated
## Tue, 14 Dec 1999 09:11:32 -0500 (EST) in response to
## Alex Ahgarin <datamanagement@email.com>
##
##


permutations <- function(n, r, v = 1:n, set = TRUE, repeats.allowed=FALSE)
{
  if(mode(n) != "numeric" || length(n) != 1 
     || n < 1 || (n %% 1) != 0) stop("bad value of n") 
  if(mode(r) != "numeric" || length(r) != 1 
     || r < 1 || (r %% 1) != 0) stop("bad value of r") 
  if(!is.atomic(v) || length(v) < n) 
    stop("v is either non-atomic or too short")
  if( (r > n) & repeats.allowed==FALSE)
    stop("r > n and repeats.allowed=FALSE")
  if(set) {
    v <- unique(sort(v))
    if (length(v) < n) stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  ## Inner workhorse
  if(repeats.allowed)
    sub <- function(n, r, v)
      {
        if(r==1) matrix(v,n,1) else
        if(n==1) matrix(v,1,r) else
        {
          inner  <-  Recall(n, r-1, v)
          cbind( rep( v, rep(nrow(inner),n)  ),
                 matrix( t(inner), ncol=ncol(inner), nrow=nrow(inner) * n ,
                        byrow=TRUE )
                )
        }
      }
  else
    sub <- function(n, r, v)
      {
        if(r==1) matrix(v,n,1) else
        if(n==1) matrix(v,1,r) else
        {
        X  <-  NULL
        for(i in 1:n)
          X  <-  rbind( X, cbind( v[i], Recall(n-1, r - 1, v[-i])))
        X
        }
      }

  sub(n, r, v[1:n])
}
