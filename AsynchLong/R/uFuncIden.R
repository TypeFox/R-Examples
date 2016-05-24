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
# xIs       : list of length nrow(data.y), the elements of which list  #
#             the rows of data.x corresponding the patient in the ith  #
#             row of data.y                                            #
#                                                                      #
#----------------------------------------------------------------------#
#                                                                      #
# Returns the value of the estimating equations.                       #
#                                                                      #
#----------------------------------------------------------------------#
uFuncIden <- function(data.y, 
                      data.x,
                      kernel, 
                      xIs,
                      yIs,
                      nPatients) {

  nCov <- ncol(data.x)

  aMat <- matrix(data = 0.0, nrow = nCov, ncol = nCov)
  bVec <- matrix(data = 0.0, nrow = nrow(data.x), ncol = nCov)
  ones <- matrix(data = 1.0, nrow = 1L, ncol = nrow(data.x))

  for( i in 1L:nPatients ) {
    ly <- yIs[[i]]$n

    if( ly < 0.5 ) next

    for( j in 1L:ly ) {

      k <- yIs[[i]]$v[j]

      lx <- xIs[[ k ]]$n

      if( lx < 0.5 ) next

      xI <- xIs[[ k ]]$v

      tx <- data.x[xI,,drop=FALSE] * kernel[[ k ]]

      bVec[xI,] <- bVec[xI,] + tx * data.y[k,3L]

      aMat <- aMat + t(tx) %*% data.x[xI,,drop=FALSE]

    }
  }

  bVec <- as.vector(ones %*% bVec)

  pars <- solve(aMat, bVec)

  return(pars)

}



