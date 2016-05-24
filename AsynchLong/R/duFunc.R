#----------------------------------------------------------------------#
# duFunc : Derivative of the estimating equation                       #
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
# kernel    : a character. One of "TD" = time dependent,               #
#             "TI" = time-invariant or "LV" = last value.              #
#             This specifies how the kernel is used                    #
#                                                                      #
# lType     : a character. One of "identity", "log", "logistic".       #
#             Specifies the form of the link function.                 #
#                                                                      #
# xIs       : list of length nrow(data.y), the elements of which list  #
#             the rows of data.x corresponding the patient in the ith  #
#             row of data.y                                            #
#                                                                      #
# yIs       : an object of class list.                                 #
#             v - which elements of y match patientID[i]               #
#             n - number of elements that match.                       #
#                                                                      #
# nPatients : an object of class numeric.                              #
#             the number of patients in dataset.                       #
#                                                                      #
#----------------------------------------------------------------------#
#                                                                      #
# Returns the value of the estimating equations.                       #
#                                                                      #
#----------------------------------------------------------------------#
duFunc <- function(pars,
                   data.y, 
                   data.x,
                   kernel, 
                   lType,
                   xIs,
                   yIs,
                   nPatients) {

  nrdx <- nrow(data.x)
  nCov <- length(pars)

  xbeta <- as.vector(data.x %*% pars)

  if( lType == "identity" ) {

    mu <- xbeta
    dmu <- data.x

  } else if( lType == "log" ) {

    exb <- exp(xbeta)
    mu <- exb
    dmu <- exb * data.x

  } else if( lType == "logistic" ) {

    exb <- exp(-xbeta)
    mu <- 1.0 / {1.0 + exb}
    dmu <- data.x * exb / {{1.0 + exb} * {1.0 + exb}}

  } else {

    stop("unsupported link function")

  }

  estEq  <- matrix(data = 0.0, nrow = nPatients, ncol = nCov)
  destEq <- matrix(data = 0.0, nrow = nCov, ncol = nCov)
  ones <- matrix(data = 1.0, nrow = 1L, ncol = nrow(data.x))

  for( i in 1L:nPatients ) {

    ly <- yIs[[i]]$n

    if( ly < 0.5 ) next

    for( j in 1L:ly ) {

      k <- yIs[[ i ]]$v[j]

      lx <- xIs[[ k ]]$n

      if( lx < 0.5 ) next

      xI <- xIs[[ k ]]$v

      tempM <- data.x[xI,,drop=FALSE] * kernel[[ k ]]

      estEq[i,] <- estEq[i,] + ones[1L,1L:lx,drop=FALSE] %*% {
                               tempM * {data.y[k,3L] - mu[xI]}}

      for( p in 1L:nCov ){
        destEq[p,] <- destEq[p,] - ones[1L,1L:lx,drop=FALSE] %*% 
                                   {dmu[xI,p] * tempM}
      }
    }
  }

  return(list(du=destEq, u=colSums(estEq), var=estEq))

}

cmp_duFunc <- compiler::cmpfun(duFunc)

