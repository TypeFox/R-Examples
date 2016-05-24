"Lcomoment.correlation" <-
function(L2) {
  if(L2$order != 2) {
    warning("L-comoment matrix argument is not of order 2")
    return()
  }

  # Following Serfling and Xiao (2006)
  #  L-correlations are the L-comoment coefficents of L-scale
  #  The diagonal of LC are the coefficients of L-variation
  LC <- Lcomoment.coefficients(L2,L2)
  return(LC)
}
