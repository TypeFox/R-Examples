calcFittedRegMat<- function( xMatAll, zMatEq, nEq, nObsEq, useMatrix, solvetol ) {
   # fitted values of regressors for IV estimations
      xMatHatEq <- list()
      for(i in 1:nEq) {
         # rows that belong to the ith equation
         rowsEq <- c( (1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i])) )
         # extract instrument matrix
         xMatAllThisEq <- xMatAll[ rowsEq, ]
         if( useMatrix ){
            xMatAllThisEq <- as( xMatAllThisEq, "dgeMatrix" )
         }
         xMatHatEq[[ i ]] <- zMatEq[[i]] %*%
            solve( crossprod( zMatEq[[i]] ),
            crossprod( zMatEq[[i]], xMatAllThisEq ), tol = solvetol )
      }
      # fitted values of all regressors
      xMatHatAll <- .stackMatList( xMatHatEq, way = "below",
         useMatrix = useMatrix )

      return( xMatHatAll )
}
