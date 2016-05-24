criticalmoisture <- maxbulkdensity <- 
function(theta, Bd, samples = NULL, graph = TRUE, ...)
{
   if (length(theta) != length(Bd))
      stop ("incompatible dimensions!")
   if (is.null(samples))
      samples <- rep(1, length(theta))
   if (length(samples) != length(theta))
      stop ("incompatible dimensions!")
   samples <- as.factor(samples)
   if (sum(tapply(theta, samples, length) < 3) > 0)
      stop ("insufficient data!")

   dat <- data.frame(theta, Bd, samples)
   lev <- levels(samples)
   nsamp <- nlevels(samples)
   lab.theta <- paste(deparse(substitute(theta)))
   lab.Bd <- paste(deparse(substitute(Bd)))
   out <- matrix(NA, nrow = 7, ncol = nsamp)
   rownames(out) <- c("Intercept", lab.theta,
      paste(lab.theta, "^2", sep = ""), "R.squared", "n",
	"critical.mois", "max.bulk")
   colnames(out) <- paste("Sample", lev, sep = " ")
   
   x <- NULL
   r.sq <- function(x, y, B)
   {
      pred <- B[1] + B[2]*x + B[3]*x^2
      sum(scale(pred, scale = FALSE)^2) / 
         sum(scale(y, scale = FALSE)^2)
   }
   
   for (j in 1:nsamp) {
      out[1:3, j] <- lm(Bd ~ theta + I(theta^2),
         data = dat[samples == lev[j], ])$coef
      out[4, j] <- r.sq(theta[samples == lev[j]],
         Bd[samples == lev[j]], out[1:3, j])
      out[6, j] <- - out[2, j] / (2 * out[3, j])
      out[7, j] <- out[1, j] + out[2, j] * out[6, j] + 
         out[3, j] * out[6, j]^2
   }
   out[5, ] <- tapply(theta, samples, length)
   
   if (graph) {
      plot(Bd ~ theta,
         pch = "",
         xlab = "Moisture",
         ylab = "Bulk density",
         main = "Soil compaction curve", ...)
      for (j in 1:nsamp) {
         points(Bd[samples == lev[j]] ~ theta[samples == lev[j]],
            col = j)
         curve(out[1, j] + out[2, j] * x + out[3, j] * x^2,
            add = TRUE, col = j)
      }
      if (nsamp > 1) {
         legend('topright', lev, lty = 1,
            col = seq(1, nsamp), cex = 0.7)
      }
   }

   out. <- list(table = out)
   class(out.) <- "criticalmoisture"
   return(out.)
}


# ------------------------------------
# print method
print.criticalmoisture <- function(x, ...)
{
    cat("\n          Critical Moisture and Maximum Bulk Density \n\n")
    print(x$table)
    invisible(x$table)
}
