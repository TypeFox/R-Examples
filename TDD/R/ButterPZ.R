ButterPZ = function(fl = NaN, nl = NaN, fh = NaN, nh = NaN, g = 1){
  poles = numeric()
  zeros = numeric()
  np = 0
  nz = 0
  Knorm = g
  if(!is.na(fl) && !is.na(nl)){
    zeros = c(zeros, rep(0, nl))
    nz = nz + nl
    poles = c(poles, 2*pi*fl*exp(1i*pi * (2*(1:nl) + nl - 1)/(2*nl)))
    np = np + nl
    # Knorm is unchanged
  }
  if(!is.na(fh) && !is.na(nh)){
    # no new zeros
    poles = c(poles, 2*pi*fh*exp(1i*pi * (2*(1:nh) + nh - 1)/(2*nh)))
    np = np + nh
    Knorm = Knorm * prod(2*pi*fh*exp(1i*pi * (2*(1:nh) + nh - 1)/(2*nh)))
  }
  return(list(Sense = 1, Knorm = Knorm, poles = signif(poles, 9), np = np, zeros = signif(zeros, 9), nz = nz)) # signifs are there to prevent small numerical errors from causing problems later
}
