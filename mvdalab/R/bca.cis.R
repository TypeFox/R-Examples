bca.cis <- function(object, conf = 0.95, 
                    type = c("coefficients", "loadings", "weights"))
{
  ncomp = 1:object$ncomp
  if ((object$val.method == "none" | object$val.method == "loo")) {
    stop("No bootstrapping was done for this model")
  }
  conf <- conf
  Bootstraps <- do.call("rbind", object$validation[names(object$validation) == type][[1]])
  Bootstraps <- as.matrix(Bootstraps[, ncomp])
  Upper <- 1 - (((1 - conf)/2))
  Lower <- 1 - Upper
  BCa.df <- llply(ncomp, function(y) {
    A <- do.call("rbind", as.list(
      by(Bootstraps[, y], list(row.names(Bootstraps)), function(x){
        Aa <- data.frame(ncomp = y, boot.mean = mean(x, na.rm = T), Skewness = skewness(x, na.rm = T))
        Aa$a <- Aa$Skewness / 6
        Aa
      }
      )))
    B <- data.frame(variable = names(Bootstraps[, 1]), Boot = Bootstraps[, y], 
                    Actual = as.vector(object[names(object) == type][[1]][, y]))
    B$Proportional.Bias <- with(B, ifelse(Boot < Actual, 1, 0))
    Proportional.Bias <- aggregate(Proportional.Bias ~ variable, B, function(x) sum(x) / object$validation$bootstraps)
    Proportional.Bias$z0 <- qnorm(Proportional.Bias$Proportional.Bias)
    Order <- as.factor(row.names(data.frame((object[names(object) == type])[[1]][, y])))
    Out <- data.frame(A, Proportional.Bias)[Order, ]
    Out
})
  BCa <- llply(ncomp, function(y) {
    Aa <- do.call("rbind", as.list(
      by(Bootstraps[, y], list(row.names(Bootstraps)), function(x){
        x
      }
    )))
    Order <- as.factor(row.names(data.frame((object[names(object) == type])[[1]][, y])))
    Aa <- Aa[Order, ]
    z0 <- BCa.df[[y]]$z0
    a <- BCa.df[[y]]$a
    a.1 <- pnorm(z0 + ((z0 - qnorm(1- ((1-conf)/2))) / (1 + a*(z0 - qnorm(1- ((1-conf)/2))))))
    a.2 <- pnorm(z0 + ((z0 + qnorm(1- ((1-conf)/2))) / (1 + a*(z0 + qnorm(1- ((1-conf)/2))))))
    Re <- do.call("rbind", llply(1:nrow(Aa), function(x) quantile(Aa[x, ], c(a.1[x], a.2[x]), na.rm = T)))
    Re <- data.frame(Re)
  names(Re) <- paste(c(Lower * 100, Upper * 100), "%", sep = "")
  Re.f <- cbind(BCa.df[[y]][, -5], Re)
  Re.f$variables <- row.names(Re.f)
  row.names(Re.f) <- NULL
  Re.f[, c(1, 9, 7:8, 5, 6, 3, 4)]
})
BCa <- sapply(ncomp, function(x) BCa[x])
BCa
}