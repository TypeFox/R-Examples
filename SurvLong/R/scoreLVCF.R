#******************************************************************************#
# Score calculation for last value carried forward method                      #
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
scoreLVCF <- function(beta, 
                      Z,  
                      X,  
                      tau, ...){

  p <- ncol(Z) - 2L

  Lmat <- matrix(data = 0.0, nrow = p, ncol = p)
  Mmatrix <- matrix(data = 0.0, nrow = p, ncol = p)
  Uvec <- numeric(length = p)

  n <- nrow(X)
  nZ <- nrow(Z)

  extractLastCov <- function(x){
    it <- which.max(x[,2L])
    x[it,]
  }

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
    # Calculate the S Function                                             #
    #----------------------------------------------------------------------#

    #----------------------------------------------------------------------#
    # Identify patients still at risk (t >= time)                          #
    #----------------------------------------------------------------------#
    ptIDs <- X[time <= X[,2L],1L]

    #----------------------------------------------------------------------#
    # Identify the covariates for this subset of patients                  #
    #----------------------------------------------------------------------#
    ZptIDs <- (Z[,1L] %in% ptIDs)

    #----------------------------------------------------------------------#
    # Keep only those covariates measured before time                      #
    #----------------------------------------------------------------------#
    ZptIDs <- ZptIDs & (Z[,2L] <= time)

    if( !any(ZptIDs) ) next

    Z2 <- Z[ZptIDs,,drop=FALSE]

    #----------------------------------------------------------------------#
    # Identify largest times and accept those covariates                   #
    #----------------------------------------------------------------------#
    cova <- by(data = Z2,
               INDICES = Z2[,1L],
               FUN = extractLastCov,
               simplify = FALSE)

    cova <- matrix(data = unlist(cova), 
                   ncol = {p+2L}, 
                   byrow = TRUE)

    IDs <- cova[,1L] == X[i,1L]
    if( sum(IDs) != 1L ) next
    cova <- cova[,c(-1L,-2L),drop=FALSE]

    prod <- exp(cova %*% beta)

    s0 <- sum(prod)

    if( (s0 > -1.5e-8) && (s0 < 1.5e-8) ) next

    s1 <- colSums(prod[,1L]*cova)

    tst <- matrix(apply(cova,1L,function(x){x %*% t(x)}),nrow=p*p)
    Zp <- rowSums(sweep(x = tst,
                        MARGIN = 2L,
                        STATS = prod[,1L],
                        FUN = "*"))
    s2 <- matrix(data = Zp, nrow = p, ncol = p)

    tmp <- cova[IDs,] - s1/s0

    Mmatrix <- Mmatrix + tmp %*% t(tmp)

    Uvec <- Uvec + tmp

    Lmat <- Lmat + (s1 %o% s1 - s2*s0)/(s0*s0)

  }

  return(list("U" = Uvec, "dUdBeta" = Lmat, "mMatrix" = Mmatrix))

}

