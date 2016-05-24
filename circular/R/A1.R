###############################################################
#                                                             #
#       R port: Claudio Agostinelli  <claudio@unive.it>       #
#                                                             #
#       Date: January, 14, 2003                               #
#       Version: 0.1-6                                        #
#                                                             #
###############################################################

A1 <- function(kappa) {
    result <- besselI(kappa, nu=1, expon.scaled = TRUE)/besselI(kappa, nu=0, expon.scaled = TRUE)
    return(result)
}
