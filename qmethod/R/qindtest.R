qindtest <- function(loa, target, nfactors) {
  # Create data frame to log the results
  qindt <- as.data.frame(matrix(NA, nrow=nfactors, ncol=2))
  qindt_log <- as.data.frame(matrix("ok", nrow=1, ncol=2))
  qindt_log[[1]] <- as.character(qindt_log[[1]])
  qindt_log[[2]] <- as.character(qindt_log[[2]])
  #qindt[1] <- as.logical(qindt[1])
  #qindt[2] <- as.logical(qindt[2])
  dimnames(qindt) <- list(paste("f",c(1:nfactors), sep=""), c("torder", "tsign"))
  dimnames(qindt_log) <- list(c("log"), c("torder_log", "tsign_log"))
  #====================================================
  # 1. Check for ORDER SWAP - REFLECTION
  #====================================================
  #Is the absolute value of diagonal coefficients larger than non-diagonal coefficients for the same factor? (if FALSE, then test is negative, OK)
  corloa <- as.data.frame(cor(loa, target))
  for (i in 1:nfactors) {
    factors <- c(1:nfactors)
    otherf <- c(factors[-i])
    #The actual test: the 'sum' is because this returns as many responses as 'other' factors, to which the given factor is compared, if one single (or more) is higher, then it is TRUE that the given factor is incorrectly positioned
    qindt[i,1] <- isTRUE(sum(abs(corloa[i,i]) < abs(corloa[otherf,i])) != 0)
  }
  if (sum(qindt[,1]) == 0) {
    qindt_log[1,1] <- "OK: No ORDER swap issues"
    warning(qindt_log[1,1])
  } else {
    #..................................................
    # Add *solution* for ORDER SWAP: reorder factors according to highest coefficients
    loa_orig <- as.matrix(loa)
    # select, for each factor, the factor that has highest correlation with the target factor
    swap <- list()
    for (j in 1:nfactors) {
      maxcor <- max(abs(corloa[,j]))
      swap <- paste(swap,which(abs(corloa[,j]) == maxcor), sep=" ") # info piece
      loa[,j] <- loa_orig[,which(abs(corloa[,j]) == maxcor)]
    }
    #..................................................
    loa <- as.matrix(loa)
    #check that no factor has been chosen twice
    corloa2 <- cor(loa)
    diag(corloa2) <- 0
    if (max(abs(corloa2)) > 0.999) {
      #!!! ehem, try format(corloa2, digits=16) and see that all correlations == 1 out of the diagonal are 0.9period instead of 1, which was causing an error!
      qindt_log[1,1] <- "ERROR in ORDER swap: at least one factor in the resample is best match for two or more factors in the target"         
      print(qindt_log[1,1])
      loa <- as.data.frame(loa_orig)
    } else {
      qindt_log[1,1] <- paste("OK - Factors reordered: ",swap[1], sep="") # (1) in the unlikely case in which the value of two or more coefficients are the same and highest, not selecting the first value would give an error---this step should actually be discarded
      if (length(swap) > 1) warning("The reordering in this iteration may be incorrect.")
      
      print(qindt_log[1,1])
      swap <- NULL
      loa_orig <- NULL
    }
  }
  #====================================================
  # 2. Check for SIGN SWAP - INVERSION
  #====================================================
  #Are all diagonal coefficients positive? (if FALSE, then test is negative, OK)
  corloa <- as.data.frame(cor(loa, target)) #remake the correlation matrix again in case that factors were reordered
  for (i in 1:nfactors) {
    qindt[i,2] <- corloa[i,i] < 0
  }
  if (sum(qindt[,2]) == 0) {
    qindt_log[1,2] <- "OK: No SIGN swap issues"
    warning(qindt_log[1,2])
    loa <- as.matrix(loa)
  } else {
    #..................................................
    # Add *solution* for SIGN SWAP: switch sign of all factor loadings
    loa_orig <- as.data.frame(loa)
    swapfactors <- which(qindt[,2] == TRUE)
    loa <- as.matrix(loa)
    loa[ , swapfactors] <- -loa[ , swapfactors]
    qindt_log[1,2] <- paste("OK: SIGN swap in factors ", paste(as.character(swapfactors), collapse=", "), sep="")
    warning(qindt_log[1,2])
    #..................................................
  }
  qindeterminacy <- list()
  qindeterminacy[[1]] <- loa
  qindeterminacy[[2]] <- qindt_log
  qindeterminacy[[3]] <- qindt
  return(qindeterminacy)
}
