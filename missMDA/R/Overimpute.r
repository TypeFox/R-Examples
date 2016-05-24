Overimpute<-function (output, plotvars = NULL) 
{
  if (length(output$call$res.Over) < 100) {
    warning("The number of imputed data sets is low to construct confidence intervals according to quantiles method. You should run MIPCA with nboot over 100.")
  }
  Overimpute_univ <- function(output, plotvars, subset = NULL) {
    X.na <- output$call$X
    if (is.null(subset)) {
      ind <- which(!is.na(X.na[, plotvars]))
    }
    else {
      ind <- subset
    }
    res <- lapply(output$call$res.Over, "[", ind, plotvars)
    res.imp <- do.call(cbind, res)
    res <- t(apply(res.imp, 1, function(x) {
      xbar <- mean(x)
      temp <- quantile(x, probs = c(0.05, 0.95))
      binf <- temp[[1]]
      bsup <- temp[[2]]
      return(c(xbar = xbar, binf = binf, bsup = bsup))
    }))
    pct <- rowMeans(is.na(X.na[ind, -plotvars]))
    res <- cbind(ind, xobs = X.na[ind, plotvars], res, pct)
    col <- cut(pct, c(-0.1, 0.2, 0.4, 0.6, 0.8, 1.1))
    levels(col) <- c("blue", "green", heat.colors(3)[c(3, 
                                                       2, 1)])
    col <- as.character(col)
    plot(x = X.na[ind, plotvars], y = res[, "xbar"], col = col, 
         xlab = "observed values", ylab = "imputed values", 
         main = colnames(X.na)[plotvars], ylim = c(min(res[,"binf"],na.rm = T), max(res[,"bsup"],na.rm = T)))
    abline(0, 1)
    segments(x0 = X.na[ind, plotvars], x1 = X.na[ind, plotvars], 
             y0 = res[, "binf"], y1 = res[, "bsup"], col = col)
    legend("topleft", legend = c("0-0.2", "0.2-0.4", "0.4-0.6", 
                                 "0.6-0.8", "0.8-1"), col = c("blue", "green", heat.colors(3)[c(3, 
                                                                                                2, 1)]), bty = "n", lty = 1, horiz = F, cex = 0.7,lwd=.4)
    return(res.over = res)
  }
  X.na <- output$call$X
  if (is.null(plotvars)) {
    Var <- 1:ncol(output$call$X)
  }
  else {
    Var <- plotvars
  }
  res <- list()
  par(mfrow = c(ceiling(sqrt(length(Var))), ceiling(length(Var)/ceiling(sqrt(length(Var))))), 
      mar = c(5, 4, 4, 2) - 1.9)
  for (i in Var) {
    subset <- which(!is.na(X.na)[, which(i == Var)])
    if (length(subset) > 1) {
      res[[as.character(i)]] <- Overimpute_univ(output = output, 
                                                plotvars = i, subset = subset)
    }
    else {
      res[[as.character(i)]] <- NULL
    }
  }
  return(res)
}
