#-------Gauss' hypergeometric function generaters-----------------------
#	my2F1.uni - generate for scalers
#	my2F1 - generate for vectors

###################################################################
# 
my2F1.uni <- function(a,b,c,z) {
  # Script to  Generate Gauss' hypergeometric function 2F1(a,b;c;z)
  # for SCALAR NEGATIVE INTEGER args a, b, and real c, z.
  n <- min(-a,-b);
  if( (n<0) || (n != round(n)) )
    stop("my2F1: Bad args: max(a, b) must be non-positive integer.");
  terms <- c(1, z*seq(-a,,-1,n)*seq(-b,,-1,n) / (seq(c,,,n)*seq(1,,,n)));
  return(sum(cumprod(terms)));
}
###################################################################

# 
my2F1 <- function(a,b,c,z) {
  # Script to  Generate Gauss' hypergeometric function 2F1(a,b;c;z)
  # for VECTOR NEGATIVE INTEGER args a, b and real c, z.
  if( (la <- length(a)) != (lb <- length(b)) ) {
    if(la==1) a <- rep(a,la<-lb);
    if(lb==1) b <- rep(b,lb<-la);
    if(la != lb)
      stop("my2F1: Bad args: a,b must be vectors of same length.");
  }
  rv <- numeric(la);
  if(length(c)<la) c  <- rep(c[1],la);
  if(length(z)<la) z  <- rep(z[1],la);
  for(i in 1:la) { rv[i] <- my2F1.uni (a[i], b[i], c[i], z[i]); }
  return(rv);
}
