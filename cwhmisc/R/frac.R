int <- function (x) {  # integer part
    as.integer(trunc(x))
}


frac <- function(x,d) {  # fractional part
  res <- abs(x-trunc(x))
  if (!missing(d)) res <- round(10^d*res)
  res
}

contfrac <- function( x, depth = 13, f=floor ) {
  ## from Mathematics 5335 Fall 2009
  ## The Euclidean Algorithm and Continued Fractions
  fac <- 10^(-3) ## to prevent Inf caused by division by nearly zero
  sgn <- sign(x)
  xi  <- x <- abs(x)
  ai <- a  <- f( xi );
  ii <- 1
  while (
         (1.0 + abs(ai-xi)*fac != 1.0) && (ii <= depth)
         ) {
    (xi <- 1.0/( xi - ai ))
    (ai <- f(xi))
    (a  <- c( a, ai) )
    ii <- ii+1
  }
  return( a*sgn )
}  # end contfrac

evalcfr <- function( cf ) {
  sgn <- sign(cf[1])
  cf  <- abs(cf)
  res <- cf[length(cf)]
  for (ii in rev(seqm(1,length(cf)-1)) ) res <-cf[ii] + 1.0/res
  return( res*sgn )
}  # end evalcfr

toCFrac <- function( x, depth=5) {
    .Machine$double.xmax^0.75
  c38 <- .Machine$double.xmax^(3/8)
  sig  <- sign(x);
  x    <- abs(x);
  minr <- 1.0/.Machine$integer.max;
  n    <- 1;
  num1 <- 1.0;      den1 <- 0.0;
  num  <- int(x);   den  <- 1.0;
  val1 <- num/den;  val2 <- 0.0;
  dif2 <- .Machine$integer.max;
  dif1 <- abs(val1-val2);
  r    <- x - num;
  while ((n < depth) & (r > c38) < (r >= minr)) {
#  while ((n < depth) & (dif1 < dif2) & (r >= minr)) {
    ## avoid division by 0 *)
   n <- n+1;
   onebyr <- 1/r + 8.*minr;  # force some rounding
   b <- int(onebyr);
   r <- onebyr - b;
   num2 <- num1;  den2 <- den1;
   num1 <- num ;  den1 <- den;
   num  <- num1*b + num2;
   den  <- den1*b + den2;
   val2 <- val1;
   val1 <- num/den;
   dif2 <- dif1;
   dif1 <- abs(val1-val2);
  } #END;  ##ELIHW*)
  return(list(num=sig*num,den=den))
} # end toCFrac

toCFrac2 <- function( x, depth=5) {
    .Machine$double.xmax^0.75
  c38 <- .Machine$double.xmax^(3/8)
  sig  <- sign(x);
  x    <- abs(x);
  minr <- 1.0/.Machine$integer.max;
  n    <- 2;
  num  <- c(1,int(x),rep(NA,depth-1))
  den  <- c(0,1,rep(NA,depth-1))
  val1 <- num[2]/den[2];  val2 <- 0.0;
  dif2 <- .Machine$integer.max;
  dif1 <- abs(val1-val2);
  r    <- x - num[2];
  while ((n <= depth) & (r > c38) < (r >= minr)) {
   n <- n+1;
   onebyr <- 1/r + 8.*minr;  # force some rounding
   b <- int(onebyr);
   r <- onebyr - b;
   num[n]  <- num[n-1]*b + num[n-2];
   den[n]  <- den[n-1]*b + den[n-2];
   val2 <- val1;
   val1 <- num[n]/den[n];
   dif2 <- dif1;
   dif1 <- abs(val1-val2);
  } #END;  ##ELIHW*)
  return(list(num=sig*num[-1],den=den[-1]))
} # end toCFrac2
