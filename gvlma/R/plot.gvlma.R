"plot.gvlma" <-
function(x, onepage = TRUE,
         ask = !onepage && prod(par("mfcol")) < 
         ncol(model.matrix(x))+4 && dev.interactive(), ...)
{
  gvlmaobj <- x
  if (!inherits(gvlmaobj, "gvlma"))
    stop("First argument must be a gvlma object.")
  X <- model.matrix(gvlmaobj)
  y <- fitted(gvlmaobj) + residuals(gvlmaobj)
  n <- nrow(X)
  yname <- all.vars(formula(gvlmaobj))[1]
  nc <- ncol(X)
  v <- gvlmaobj$GlobalTest$timeseq
  fits <- fitted(gvlmaobj)
  sighat <- sqrt(sum(resid(gvlmaobj)^2)/n)
  stdresids <- residuals(gvlmaobj)/sighat
  if (onepage)
    {
      op <- par(mfrow = c((6 + nc-1)/2, 2), col = 1)
      on.exit(par(op))
    }
  if (ask)
    {
      op <- par(ask = TRUE)
      on.exit(par(op))
    }
  for(j in 2:nc) {
    plot(X[, j], y, xlab = labels(X)[[2]][j], main = 
         "Plot of Response Variable versus Predictor Variable",
         ylab = yname, lty = 1, type = "p")
  }
  plot(v, y, main = 
       "Plot of Response Variable versus Time Sequence", xlab
       = "Time Sequence for Directional Stat 4", ylab = yname, lty = 1, type = 
       "p")
  plot(fits, stdresids, main = 
       "Plot of the Standardized Residuals versus the Fitted Values",
       xlab = "Fitted Values", ylab = "Standardized Residuals"
       )
  hist(stdresids, col="blue", main = 
       "Histogram of the Standardized Residuals", xlab = 
       "Standardized Residuals", ylab = "Frequency")
  qqnorm(stdresids, main = 
         "Normal Probability Plot of the Standardized Residuals (with line)"
         )
  qqline(stdresids)
  plot(v, stdresids, main = 
       "Plot of the Standardized Residuals versus Time Sequence",
       xlab = "Time Sequence for Directional Stat 4",
       ylab = "Standardized Residuals" 
       )  
}

