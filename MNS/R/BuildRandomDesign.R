BuildRandomDesign <-
function(fixedDes, blups, N){
  # build random design matrix based on fixed effects design and blups
  #
  # This is as described in Bondell et al. (2010)
  #
  # INPUT:
  #      - fixedDes: design matrix for fixed effects. The design for the random effects will be the same expect multipled by corresponding BLUPs
  #      - blups: vector of blups
  #      - N: number of subjects (so we know which data corresponds to whom)
  #
  #
  
  p = ncol(fixedDes)
  nobs = nrow(fixedDes)/N # number of observations per subject
  nsub = N
  
  for (sub in 1:N){
    ID = seq(1+(sub-1)*p,(sub*p))
    ID2 = seq(1+(sub-1)*nobs,(sub*nobs))
    fixedDes[ID2, ] = fixedDes[ID2, ] %*% diag(blups[ID])
  }
  
  return(fixedDes)
}
