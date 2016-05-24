#******************************************************************************#
# Score calculation for half kernel method                                     #
#******************************************************************************#
#                                                                              #
# Inputs                                                                       #
#                                                                              #
#  beta           an object of class numeric.                                  #
#                 parameter estimate(s)                                        #
#                                                                              #
#  Z              an object of class data.frame.                               #
#                 The structure of the data.frame must be                      #
#                 \{patient ID, time of measurement, measurement(s)\}.         #
#                 Patient IDs must be of class integer or be able to be        #
#                 coerced to class integer without loss of information.        #
#                 Missing values must be indicated as NA.                      #
#                                                                              #
#  X              an object of class data.frame.                               #
#                 The structure of the data.frame must be                      #
#                 \{patient ID, event time, event indicator\}.                 #
#                 Patient IDs must be of class integer or be able to be        #
#                 coerced to class integer without loss of information.        #
#                 Missing values must be indicated as NA.                      #
#                                                                              #
#  tau            an object of class numeric.                                  #
#                 The desired time point.                                      #
#                                                                              #
#  h              an object of class numeric.                                  #
#                 kernel bandwidth                                             #
#                                                                              #
#  kType          an object of class character indicating the type of          #
#                 smoothing kernel to use in the estimating equation.          #
#                 Must be one of \{"epan", "uniform", "gauss"\}, where         #
#                 "epan" is the Epanechnikov kernel and "gauss" is the         #
#                 Gaussian kernel.                                             #
#                                                                              #
#  Outputs                                                                     #
#                                                                              #
#  Returns a list                                                              #
#                                                                              #
#  U              an object of class numeric.                                  #
#                 Score function(s)                                            #
#                                                                              #
#  dUdBeta        an object of class numeric matrix                            #
#                 Derivative of the Score function                             #
#                                                                              #
#  mMatrix        an object of class numeric matrix                            #
#                 Sigma                                                        # 
#                                                                              #
#******************************************************************************#
scoreHalf <- function(beta, 
                      Z, 
                      X,  
                      tau,  
                      h,
                      kType, ...){

  p <- ncol(Z) - 2L

  Zp <- data.matrix(Z[,3L:{p+2L},drop=FALSE])
  Zid <- Z[,1L]
  Ztime <- Z[,2L]
  rm(Z)

  dUdBeta <- matrix(data = 0.0, nrow = p, ncol = p)
  Mmatrix <- matrix(data = 0.0, nrow = p, ncol = p)
  Uvec <- numeric(length = p)

  n <- nrow(X)

  outerOnce <- matrix(data = apply(Zp,1L,function(x){x %*% t(x)}),
                      nrow = nrow(Zp),
                      ncol = p*p,
                      byrow = TRUE)

  for( i in 1L:n ) {

    #----------------------------------------------------------------------#
    # If the time is censored, do not include in summation                 #
    #----------------------------------------------------------------------#
    if( X[i,3L] < 0.5 ) next

    #----------------------------------------------------------------------#
    # If the time is greater than the integration limit, do not include in #
    # the summation                                                        #
    #----------------------------------------------------------------------#
    time <- X[i,2L]

    if( time > tau ) next

    #----------------------------------------------------------------------#
    # Keep only those covariates with measurement time less equal time     #
    #----------------------------------------------------------------------#
    use <- (Ztime <= time) & Zid == X[i,1L]

    if( !any(use) ) next

    #----------------------------------------------------------------------#
    # Calculate the S Function                                             #
    #----------------------------------------------------------------------#
    kern <- local_kernel(t = {time - Ztime}, 
                         h = h,  
                         kType = kType )

    #----------------------------------------------------------------------#
    # Identify patients still at risk (t >= time)                          #
    #----------------------------------------------------------------------#
    ptIDs <- X[time <= X[,2L],1L]

    #----------------------------------------------------------------------#
    # Identify the covariates for this subset of patients                  #
    #----------------------------------------------------------------------#
    ZptIDs <- (Zid %in% ptIDs)

    #----------------------------------------------------------------------#
    # Keep only those covariates measured before time                      #
    #----------------------------------------------------------------------#
    ZptIDs <- ZptIDs & (Ztime <= time)

    cova <- Zp[ZptIDs,,drop = FALSE]

    prod <- kern[ZptIDs]*exp(cova %*% beta)

    s0 <- sum(prod)

    if( (s0 > -1.5e-8) && (s0 < 1.5e-8) ) next

    s1 <- colSums(prod[,1L]*cova)

    s2 <- matrix(data = colSums(outerOnce[ZptIDs,,drop=FALSE]*prod[,1L]),
                 nrow = p,
                 ncol = p)

    #----------------------------------------------------------------------#
    # Calculate U and dUdBeta                                              #
    #----------------------------------------------------------------------#
    ZmZ <- sweep(x = Zp[use,,drop = FALSE], 
                 MARGIN = 2L,
                 STATS = s1/s0,
                 FUN = "-")

    tmp <- colSums( kern[use] * ZmZ )

    Mmatrix <- Mmatrix + tmp %*% t(tmp)

    Uvec <- Uvec + tmp

    dUdBeta <- dUdBeta + sum(kern[use])*(s1 %*% t(s1) - s2*s0)/(s0*s0)

  }

  return(list("U" = Uvec, "dUdBeta" = dUdBeta, "mMatrix" = Mmatrix))

}

