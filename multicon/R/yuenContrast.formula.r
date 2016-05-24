yuenContrast.formula <- function(formula, data = NULL, wgt = c(1, -1), tr = .2, alpha = .05, EQVAR = FALSE, alternative = "unequal", ...) {
  if(tr > .25) print("Warning: with tr > .25 type I error control might be poor.")
  z <- model.frame(formula, data = data)
  Z <- na.omit(z)
  GrpNs <- tapply(Z[,1], Z[,2], length)
  if(length(GrpNs)!=length(wgt)) {stop("Each level of the IV must have a specified contrast weight")}
  GrpMs <- tapply(Z[,1], Z[,2], mean, tr=tr)
  GrpVars <- tapply(Z[,1], Z[,2], winvar, tr=tr)
  GrpHs <- GrpNs - floor(tr*GrpNs)*2
  GrpDs <- ((GrpNs - 1)*GrpVars) / (GrpHs*(GrpHs-1))
  if(EQVAR==FALSE) {
    DF <- (sum(GrpDs)^2) / sum(GrpDs^2 / (GrpHs-1))
    Num <- sum(GrpMs*wgt)
    Den <- sqrt(sum(GrpVars*wgt^2 * ((GrpNs - 1)/(GrpHs*(GrpHs-1)))))
    teststat <- Num / Den
  }
  if(EQVAR==TRUE) {
    DF <- sum(GrpHs - 1)
    PoolVar <- sum((GrpNs - 1)*GrpVars) / sum(GrpNs-1)
    Num <- sum(GrpMs*wgt)
    Den <- sqrt(sum(wgt^2*((GrpNs - 1)/(GrpHs*(GrpHs-1))))*PoolVar) # Keep checking this line
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
  Ms <- data.frame(GrpNs, GrpMs, wgt)
  colnames(Ms) <- c("N", "M", "wgt")
  rownames(Ms) <- names(GrpMs)
  test <- data.frame("stat"=teststat, "df"=DF, "crit"=crit, "p"=p, row.names="results")
  out <- list(Ms, test)
  names(out) <- c("Ms", "test")
  return(out)
}
