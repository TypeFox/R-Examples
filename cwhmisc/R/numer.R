scm <- function( m, n ) {n*m/gcd(m,n)}  ## Smallest Common Multiple

EulerPhi <- function(n) {sum(unlist(lapply(1:n,function(x) gcd(x,n)==1)))}

gcd <- function( a, b ) {
  Euclid(a,b)[3]
}  ## end EulerPhi

Euclid <- function( a, b ) {
  ## http://en.wikipedia.org/wiki/Extended_Euclid%27s_algorithm
  S <- 0;    OldS <- 1
  T <- 1;    OldT <- 0
  R <- min(abs(c(a,b)));  OldR <- max(abs(c(a,b)))
  while (R != 0) {
      q <- OldR %/% R
      prov <- R;  R <- OldR - q*R; OldR <- prov
      prov <- S;  S <- OldS - q*S; OldS <- prov 
      prov <- T;  T <- OldT - q*T; OldT <- prov 
  }
  return (c(OldS, OldT,OldR))
}  ## end Euclid

Inv <- function(a, n) {
  T <- 0;    NewT <- 1
  R <- abs(n);    NewR <- abs(a)
  while (NewR != 0) {
    q <- R %/% NewR
    prov <- NewT;  NewT <- T - q*NewT; T <- prov 
    prov <- NewR;  NewR <- R - q*NewR; R <- prov
  }
  if ( R > 1)  return(NA) else {
    if( T <  0 ) T <- T + n
    return (T)
  }
}  # end Inv

modexp <- function(a, b, n)  {  ## a^b mod n  using repeated squaring
  ## from http://mvngu.wordpress.com/2008/08/01/parigp-programming-for-basic-cryptography/
    bin <- intToBase( b )
    d <- a %% n;
    for (ii in seqm(2,nchar(bin))) {
        d <- (d*d) %% n
        if (substr(bin,ii,ii) == "1")  d <- (d*a) %% n
    }
    return(d %% n);
}  ## end modexp
