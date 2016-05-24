################################################################################
#  Confidence intervals of regression parameters' estimates                    #
################################################################################
#                                                                              #
#  Computes confidence intervals of regression parameters' estimates           #
#  for objects of class 'parfm'                                                #
#                                                                              #
#  Its parameters are                                                          #
#   - x         : the fitted model, object of class 'parfm'                    #
#   - level     : the coverage probability of the interval                     #
#   - digits    : number of significant digits                                 #
#                                                                              #
#                                                                              #
#   Date: January, 16, 2012                                                    #
#   Last modification on: February, 29, 2012                                  #
################################################################################

ci.parfm <- function(x,
                     level=.05,
                     digits=3) {
  beta <- which (!is.na(x[, "p-val"]))
  
  q <- qnorm(1 - level / 2)
  
  res <- exp(x[beta, "ESTIMATE"] + outer(x[beta, "SE"], c(-1, 1) * q))

  colnames(res) <- c("low", "up")
  rownames(res) <- names(beta)
  
  return(round(res, digits))
}
  