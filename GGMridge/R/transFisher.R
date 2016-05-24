#----------------------------------------------------------------------#
# Fisher's Z-transformation of (partial) correlation.                  #
#----------------------------------------------------------------------#
# Inputs :                                                             #
#                                                                      #
#  x  A vector having entries between -1 and 1                         #
#                                                                      #
# Outputs :                                                            #
#                                                                      #
#  Fisher's Z-transformed values.                                      #
#                                                                      #
#----------------------------------------------------------------------#
transFisher <-  function(x) {

  eps <- 1e-7

  x <- pmin(x,  1 - eps)

  x <- pmax(x, -1 + eps)

  fish <- 0.5 * log( {1.0 + x} / {1.0 - x} )

  return( fish )

}
