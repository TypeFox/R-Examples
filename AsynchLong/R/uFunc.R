#----------------------------------------------------------------------#
# uFunc : Estimating equation                                          #
#----------------------------------------------------------------------#
#                                                                      #
# pars      : Current parameter estimates                              #
#                                                                      #
# data.y    : Matrix of responses. It is assumed that the first column #
#             contains integer patient IDs, the second column contains #
#             the time of measurement, and the third column contains   #
#             the value of the measurement.                            #
#                                                                      #
# data.x    : Matrix of covariates. The columns contain only the values#
#             of the covariates.                                       #
#                                                                      #
# kernel    : a list, ith element containing a matrix of the yIs       #
#             distances                                                #
#                                                                      #
# lType     : a character. One of "identity", "log", "logistic".       #
#             Specifies the form of the link function.                 #
#                                                                      #
# xIs       : list of length nrow(data.y), the elements of which list  #
#             the rows of data.x corresponding the patient in the ith  #
#             row of data.y                                            #
#                                                                      #
#----------------------------------------------------------------------#
#                                                                      #
# Returns the value of the estimating equations.                       #
#                                                                      #
#----------------------------------------------------------------------#
uFunc <- function(pars,
                  data.y, 
                  data.x,
		  kernel, 
                  lType,
                  xIs,
                  yIs,
                  nPatients) {

  xbeta <- data.x %*% pars

  if( lType == "logistic" ) {
    mu <- 1.0/(1.0 + exp(-xbeta))
  } else {
    mu <- exp(xbeta)
  }

  y3 <- data.y[, 3L]

  zeroVec <- numeric(length=length(pars))

  tempFunc <- function(x){

    if( xIs[[ x ]]$n < 0.5 ) return( zeroVec )

    xI <- xIs[[ x ]]$v

    res <-  colSums( data.x[xI,,drop=FALSE] * kernel[[ x ]] * 
                     {y3[x] - mu[xI]} )

    return(res)
  }

  tempFunc2 <- function(x){

    if( yIs[[x]]$n < 0.5 ) return( zeroVec )

    temp <- rowSums(sapply(yIs[[x]]$v, tempFunc))

    return(temp)

  }

  estEq <- rowSums(sapply(1L:nPatients, tempFunc2))

  return(estEq)

}

cmp_uFunc <- compiler::cmpfun(uFunc)


