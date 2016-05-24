tContrast.default <- function(IV, DV, wgt = c(1, -1), alpha = .05, EQVAR = FALSE, alternative = "unequal", ...) {
  if(length(unique(IV))!=length(wgt)) {stop("Each level of the IV must have a specified contrast weight")}
  z <- data.frame(DV, IV)
  Z <- subset(z, complete.cases(z))
  GrpNs <- tapply(Z[,1], Z[,2], length)
  GrpMs <- tapply(Z[,1], Z[,2], mean)
  GrpVars <- tapply(Z[,1], Z[,2], var)
  if(EQVAR==FALSE) {
    GrpQs <- GrpVars / GrpNs
    DF <- sum(GrpQs)^2 / sum(GrpQs^2 / (GrpNs-1))
    Num <- sum(GrpMs*wgt)
    Den <- sqrt(sum((GrpVars*wgt^2)/GrpNs))
    teststat <- Num / Den
  }
  if(EQVAR==TRUE) {
    DF <- sum(GrpNs-1)
    PoolVar <- sum((GrpNs - 1)*GrpVars) / sum(GrpNs-1)
    Num <- sum(GrpMs*wgt)
    Den <- sqrt(sum(wgt^2/GrpNs)*PoolVar)
    teststat <- Num / Den
  }
  crit <- ifelse(alternative=="unequal", qt(1-alpha/2,DF), qt(1-alpha,DF))
  if(alternative=="unequal") {
    p <- 2*(1-pt(abs(teststat),DF))
  }
  if(alternative=="greater") {
    p <- ifelse(teststat > 0, 1-pt(abs(teststat),DF), ifelse(teststat < 0, pt(abs(teststat),DF), .5))
  }
  if(alternative=="less") {
    p <- ifelse(teststat > 0, pt(abs(teststat),DF), ifelse(teststat < 0, 1-pt(abs(teststat),DF), .5))
  }
  r.contrast <- sqrt(teststat^2 / (teststat^2 + sum(GrpNs-1)))
  Ms <- data.frame(GrpNs, GrpMs, wgt)
  colnames(Ms) <- c("N", "M", "wgt")
  rownames(Ms) <- names(GrpMs)
  test <- data.frame("stat"=teststat, "df"=DF, "crit"=crit, "p"=p, "r_contrast"=r.contrast, row.names="results")
  out <- list(Ms, test)
  names(out) <- c("Ms", "test")
  return(out)
}
