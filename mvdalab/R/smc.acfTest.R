smc.acfTest <- function(object, ncomp = object$ncomp) {
  if(length(ncomp) > 1) {
    stop("Sorry...assessment is one LV at-a-time")
  }
X.error.smc <- smc(object, ncomp)$error
  acf.ci <- function(COLUMN, ci) {
    clim0 <- qnorm((1 + ci) / 2)/sqrt(length(COLUMN))
    clim <- clim0 * sqrt(cumsum(c(1, 2 * acf(COLUMN, plot = F)$acf[-1, 1, 1]^2)))
    ACF.LB <- max(-clim)
    ACF.UB <- min(clim)
  return(c(ACF.LB, ACF.UB))
  }
  ACF <- Significant <- NULL
  ACFs <- llply(colnames(X.error.smc), function(this.col) {
    Sig.ACFs <- data.frame(variable = this.col, ACF = acf(X.error.smc[, this.col], plot = F)$acf[-1, 1, 1], 
        Significant = ifelse(acf(X.error.smc[, this.col], plot = F)$acf[-1, 1, 1] > acf.ci(X.error.smc[, this.col], .95)[1] &
        acf(X.error.smc[, this.col], plot = F)$acf[-1, 1, 1] < acf.ci(X.error.smc[, this.col], .95)[2], 
              "Not.Significant", "Significant"))
    Sig.ACFs.F <- subset(Sig.ACFs, Significant == "Significant" & ACF > 0)
    Sig.ACFs.F2 <- Sig.ACFs.F[diff(c(0, as.numeric(rownames(Sig.ACFs.F)))) == 1, ][, 1]
      if(length(Sig.ACFs.F2) == 0) {
        K <- 0
      } else {
        K <- 1:length(Sig.ACFs.F2)
      } 
      Sig.ACFs[K, ]
     })
ACFs
}