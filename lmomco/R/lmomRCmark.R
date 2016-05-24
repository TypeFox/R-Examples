"lmomRCmark" <-
function(x, rcmark=NULL, r=1, sort=TRUE) {
  n <- length(x);
  if(is.null(rcmark)) rcmark <- rep(0,n);
  if(n != length(rcmark))
     stop("sample size does not match size of censoring indicator");

  if(sort) {
    ix <- sort(x, index.return=TRUE)$ix;
    x <- x[ix]; rcmark <- rcmark[ix];
  }
  xn <- x[n];

  "fnG" <- function(t) {
     if(t >= xn) return(0);
     xs <- x[x <= t];   cs <- rcmark[x <= t];
     nxs <- length(xs); js <- 1:nxs;
     G <- ( (n - js) / (n - js + 1) )^(1-cs);
     return(cumprod(G)[nxs]);
  }

  "wj" <- function(r, j) {
     ssum <- 0;
     for(k in 0:(r-1)) {
        G.right <- fnG(x[j]);
        ifelse((j-1) == 0, G.left <- 1, G.left <- fnG(x[j-1]));
        beta.left  <- pbeta(1 - G.right, shape1=r-k, shape2=k+1);
        beta.right <- pbeta(1 - G.left,  shape1=r-k, shape2=k+1);
        ssum <- ssum + (-1)^k*choose(r-1,k)*(beta.left - beta.right);
     }
     return(ssum/r);
   }

  weight.values <- sapply(1:n, wj, r=r);
  return(sum(weight.values*x));
}




