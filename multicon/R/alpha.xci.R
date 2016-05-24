alpha.xci <-
function(x, k, n, CI=.95) {
  pcrit <- (1-CI)/2
  df1 <- n*(k-1)
  df2 <- n
  fcrit.L <- qf(pcrit,df1,df2,lower.tail=TRUE)
  fcrit.U <- qf(pcrit,df1,df2,lower.tail=FALSE)
  LL <- 1-((1-x)/fcrit.L)
  UL <- 1-((1-x)/fcrit.U)
  out <- cbind(LL, UL)
  colnames(out) <- c("Lower Limit", "Upper Limit")
  return(out)
}
