#flags Q sorts automatically according to the given loadings matrix
qflag <- function(loa=loa, nstat) {
  # calculate number of Q sorts and number of statements
  nqsorts <- nrow(loa)
  #FLAGGING CRITERIA: 
  # -- 1) qsorts which factor loading is higher than the threshold for pval >0.95, and 
  # -- 2) qsorts which square loading is higher than the sum of square loadings of the same q-sort in all other factors
  thold.05 <- 1.96/sqrt(nstat)
  loa_sq <- loa^2
  flagged <- data.frame(cbind(1:nqsorts))
  flagged[,1] <- as.logical(flagged[,1])
  f <- 1
  while (f <= ncol(loa)) {
    n <- 1
    while (n <= nqsorts) {
      flagged[n,f] <- loa_sq[n,f] > (rowSums(loa_sq)[[n]]-loa_sq[n,f]) & abs(loa[n,f]) > thold.05
      n <- n+1
    }
    f <- f+1
  }
  flagged <- as.matrix(flagged)
  colnames(flagged) <- paste("flag_f",1:ncol(loa), sep="")
  row.names(flagged) <- row.names(loa)
  return(flagged)
}