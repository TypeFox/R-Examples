WeightsRootsFejer <- function(n){
  
  ## Adpoted from FEJER2_RULE_COMPUTE by Joerg Waldvogel
  # Reference:
  #    Joerg Waldvogel,
  #    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
  #    BIT Numerical Mathematics
  #    Volume 43, Number 1, pages 1-18, 2003.
  
  x <- rep(NA,n)
  w <- rep(NA,n)
  
  if ( n == 1 ){
    x[1] <- 0.0;
    w[1] <- 2.0;    
  } else {
    
    N = seq(from=1,to=n,by=2);
    L = length ( N );
    m = n + 1 - L;
  
    v0 = c(2/N/(N-2),1/N[L],rep(0,m));
    v2 = -v0[1:(n+1)] - v0[seq(from=n+2,to=2,by=-1)];
  
    w = Re(fft(v2,inverse=TRUE)/length(v2));
    x = cos ( pi * ( seq(from=n,by=-1,to=1) / ( n + 1 ) ));

    w = w[2:(n+1)];
  }
    
  res <- list(Xn=x,Wn=w)
  return(res)
}