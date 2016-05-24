# Define R function to compute Nonparametric Multivariate Change Point Model #
NPMVCP <- function(X) {
  X <- as.matrix(X)

  # dimension of each data vector is p = number of columns in X #
  p <- ncol(X)

  # set quarantine constant based on p #
  if (p == 2) { c <- 9 } else { c <- 15 }

  # N = total number of observation vectors present in X #
  N <- nrow(X)

  # set monitoring start value #
  monitoring.start <- max(p + 10, 2*c + 3)

  # initalize matrix to hold R_n(X^(i)) #
  Rvecsn <- matrix(0, nrow=N, ncol=p)
  # initialize matrix to hold Di,n #
  Din <- matrix(0, nrow=N, ncol=p)
  # initialize vector to hold maxRkn #
  maxRkn <- rep(0, N - (monitoring.start-1))
  # initialize vector to hold location of maximum #
  maxRknloc <- rep(NA, N - (monitoring.start-1))

  # compute Rkn test statistic for each n >= monitoring.start sequentially #
  for (n in 1:N) {
    # update centered rank vectors #
    for (i in 1:n){
      diff <- X[i,] - X[n,]
      denom <- sqrt(sum(diff^2))
      Din[i,] <- rep(0, p)
      if (denom != 0) { Din[i,] <- diff/denom }
      if (i < n){ Rvecsn[i,] <- Rvecsn[i,] + Din[i,] }
    }
    # Note: Dn,i = - Di,n #
    # Also, if n = 1, R_1(X_1) = (0, dotsc, 0)^T #
    if (n > 1){
      Rvecsn[n,] <- apply(-Din[1:n,], 2, sum)
    }

    if (n >= monitoring.start) {
      SIGMAhatn <- (1/(n-1))*(t(Rvecsn[1:n,]) %*% Rvecsn[1:n,])
      SIGMAhatninv <- solve(SIGMAhatn)
      for (k in 1:(n-1)){
        if (k == 1){
          Rbarnk <- Rvecsn[1,]
        } else {
          Rbarnk <- ((k-1)/k)*Rbarnkmin1 + (1/k)*Rvecsn[k,]
        }
        SIGMAhatRkninv <- (n*k/(n-k))*SIGMAhatninv
        Rkn <- t(Rbarnk) %*% SIGMAhatRkninv %*% Rbarnk
        Rbarnkmin1 <- Rbarnk
        if ((Rkn > maxRkn[n-(monitoring.start-1)]) & (c < k) & (k < (n-c))) {
          maxRkn[n-(monitoring.start-1)] <- Rkn
          maxRknloc[n-(monitoring.start-1)] <- k
        }
      }
    }
  }

  # create data frame for output #
  output <- data.frame(X)
  output$Rmax <- c(rep(NA, monitoring.start-1), maxRkn)
  output$tauhat <- c(rep(NA, monitoring.start-1), maxRknloc)
  return(output)
}
