#----------------------------------------------------------------------#
# Given a geno or nongeno object, calculate the kernel.                #
#----------------------------------------------------------------------#
calcKernObj <- function(obj) {

  wgts <- obj@weights

  if( obj@kernel == "ibs" ) {

    #--------------------------------------------------------------#
    # Calculate IBS kernel.                                        #
    #--------------------------------------------------------------#
    wgts <- as.vector(obj@weights)

    deno <- sum(wgts)*2.0

    wgts <- matrix(data = sqrt(wgts), 
                   nrow = nrow(obj@mat), 
                   ncol = ncol(obj@mat), 
                   byrow = TRUE)

    ind0 <- 1.0*(obj@mat < 0.5) * wgts
    ind1 <- 1.0*({obj@mat > 0.5} & {obj@mat < 1.5}) * wgts
    ind2 <- 1.0*({obj@mat > 1.5} & {obj@mat < 2.5}) * wgts

    ind <- tcrossprod({2.0*ind0 + ind1}, ind0) + 
           tcrossprod({2.0*ind1 + ind0 + ind2}, ind1) +
           tcrossprod({2.0*ind2 + ind1}, ind2)

    res <- ind/deno

  } else if( obj@kernel == "linear" ) {

    #--------------------------------------------------------------#
    # Calculate linear kernel.                                     #
    #--------------------------------------------------------------#
    if( ncol(wgts) == 1L ) {
      wgts <- diag(wgts[,1L], nrow = nrow(wgts), ncol = nrow(wgts))
    }

    res <- obj@mat %*% wgts %*% t(obj@mat)

  } else if( obj@kernel == "quadratic" ) {

    #--------------------------------------------------------------#
    # Calculate quadratic kernel.                                  #
    #--------------------------------------------------------------#
    if( ncol(wgts) == 1L ) {
      wgts <- diag(wgts[,1L], nrow = nrow(wgts), ncol = nrow(wgts))
    }

    res <- ( 1.0 + obj@mat %*% wgts %*% t(obj@mat))^2
 
  } else if( obj@kernel == "interactive" ) {

    #--------------------------------------------------------------#
    # Calculate interactive kernel.                                #
    #--------------------------------------------------------------#
    res <- 1.0 + obj@mat %*% t(obj@mat)

    k <- 2L
    while( k <= ncol(obj@mat) ) {
      temp <- obj@mat[,1L:{k-1L},drop=FALSE] * obj@mat[,k]
      res <- res + temp %*% t(temp)
      k <- k + 1L
    }

  } else {

    stop("kernel must be one of {ibs, linear, quadratic, interactive}.",
         call. = FALSE)

  }

  return(res)

}

