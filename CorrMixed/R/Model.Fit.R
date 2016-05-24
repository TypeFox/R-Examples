Model.Fit <- function(Model.1, Model.2){
 
if ((Model.1$Model=="Model 1, Random intercept" & Model.2$Model=="Model 2, Random intercept + serial corr (Gaussian)")|
  (Model.2$Model=="Model 1, Random intercept" & Model.1$Model=="Model 2, Random intercept + serial corr (Gaussian)")){
  
  G2 <- 2 * abs(Model.1$LogLik - Model.2$LogLik)
  
  p <- pchisq(c(G2), df=2, lower.tail=FALSE)
}


if ((Model.1$Model=="Model 2, Random intercept + serial corr (Gaussian)" & Model.2$Model=="Model 3, Random intercept, slope + serial corr (Gaussian)")|
      (Model.2$Model=="Model 3, Random intercept, slope + serial corr (Gaussian)" & Model.1$Model=="Model 2, Random intercept + serial corr (Gaussian)")){
  
  G2 <- 2 * abs(Model.1$LogLik - Model.2$LogLik)
  
  p <- .5 * pchisq(c(G2), df=1, lower.tail=FALSE) +
    .5 * pchisq(c(G2), df=2, lower.tail=FALSE)
  }

cat("Compare Models:\n")
cat("=================\n")
cat("G2: ", G2, sep="")
cat("\np-value: ", p, sep="")

}

