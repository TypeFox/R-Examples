mag2flux=function( mag, zero_pt=21.10, ABwave=F) {

  if(!missing(ABwave))
    return(10^(-0.4*(mag + 2.406 + 5*log10(ABwave))))
  else 
    return(10^(-0.4*( mag + zero_pt)))
}
