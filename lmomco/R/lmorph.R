"lmorph" <-
function(lmom) {
  # if any of the trimming characteristics are non zero then prevent
  # conversion to the named attribute L-moment object
  if((length(lmom$trim)     == 1 && lmom$trim     > 0 ) ||
     (length(lmom$leftrim)  == 1 && lmom$leftrim  > 0 ) ||
     (length(lmom$rightrim) == 1 && lmom$rightrim > 0 ) ) {
    warning("L-moment argument appears to contain non zero trimming. I can not morph to alternative L-moment object")
    return(lmom)
  }


  if(length(lmom$L1) == 1) {
    # object is like the lmom.ub object--thus, trimming is 0
    L <- c(lmom$L1,lmom$L2,lmom$L3,lmom$L4,lmom$L5)
    R <- c(NA,lmom$LCV,lmom$TAU3,lmom$TAU4,lmom$TAU5)
    z <- list(lambdas = L, ratios = R,
              trim=0, leftrim=NULL, rightrim=NULL, source="lmorph")
  }
  else {
    # note that higher L-moments than order 5 are thrown away
    z <- list(L1   = lmom$lambdas[1],
              L2   = lmom$lambdas[2],
              TAU3 = lmom$ratios[3],
              TAU4 = lmom$ratios[4],
              TAU5 = lmom$ratios[5],
              LCV  = lmom$ratios[2],
              L3   = lmom$lambdas[3],
              L4   = lmom$lambdas[4],
              L5   = lmom$lambdas[5],
              source = "lmorph")
  }
  return(z)
}
