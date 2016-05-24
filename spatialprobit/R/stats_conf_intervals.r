# PURPOSE: Computes an hperc-percent credible interval for a vector of MCMC draws
# --------------------------------------------------------------------
# Usage: bounds = cr_interval(draws,hperc);
# where draws = an ndraw by nvar matrix
#       hperc = 0 to 1 value for hperc percentage point
# --------------------------------------------------------------------
# RETURNS:
#         bounds = a 1 x 2 vector with 
#         bounds(1,1) = 1-hperc percentage point
#         bounds(1,2) = hperc percentage point
#          e.g. if hperc = 0.95
#          bounds(1,1) = 0.05 point for 1st vector in the matrix
#          bounds(1,2) = 0.95 point  for 1st vector in the matrix
#          bounds(2,1) = 0.05 point for 2nd vector in the matrix
#          bounds(2,2) = 0.05 point for 2nd vector in the matrix
#          ...
# --------------------------------------------------------------------
# Written by J.P. LeSage
# This function takes a vector of MCMC draws and calculates
# an hperc-percent credible interval
#
# Adpted to R by
# Miguel Godinho de Matos
# Dept. Engineering & Public Policy
# Carnegie Mellon University
mcmc_conf_interval <- function(adraw, hperc){
  ndraw  <- nrow( adraw )
  nvars  <- ncol( adraw )
  lIdx   <- round( (0.5 - hperc/2 ) * ndraw, 0 )   #lower index in the vector
  uIdx   <- round( (0.5 + hperc/2 ) * ndraw, 0 )   #upper index in the vector
  bounds <- matrix(data = 0 , nrow=nvars, ncol=2 ) #solution container
  for( i in 1:nvars ){
    tmp        <- sort( adraw[,i] )
    bounds[i,] <- c( tmp[lIdx] , tmp[uIdx] ) 
  }
  return( bounds )
}

