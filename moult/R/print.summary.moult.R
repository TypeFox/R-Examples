print.summary.moult <- function(x, digits = max(3, getOption("digits") - 3), ...)
 {  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 
        0.85)), "", sep = "\n")
    if (!x$converged) {
      cat("model did not converge\n")
    }
    else {
    
      n.d <- length(x$coefficients$duration)
      n.m <- length(x$coefficients$mean)
      n.s <- length(x$coefficients$sd)
      n.par <- n.d + n.m + n.s

      coefmat <- data.frame(unlist(x$coefficients), unlist(x$standard.errors))
      names(coefmat) <- c("Estimate", "Std. Error") 
  #    rownames(coefmat) <- c(names(x$coefficients$duration), names(x$coefficients$mean),
  #                           names(x$coefficients$sd))

      cat("\nDuration coefficients:\n")

      d.mat <- coefmat[(1:n.d), ]
      rownames(d.mat) <- names(x$coefficients$duration)
      printCoefmat(d.mat, digits = digits, signif.legend = FALSE)

      cat("\nMean start date coefficients:\n")

      m.mat <- coefmat[(n.d + 1):(n.d + n.m), ]
      rownames(m.mat) <- names(x$coefficients$mean)
      printCoefmat(m.mat, digits = digits, signif.legend = FALSE)

      cat("\nCoefficients for standard deviation in start date:\n")

      sd.mat <- coefmat[(n.d + n.m + 1):(n.par), ]
      rownames(sd.mat) <- names(x$coefficients$sd)
      printCoefmat(sd.mat, digits = digits, signif.legend = FALSE)

      cat("\nLog-likelihood:", formatC(x$loglik, digits = digits), 
            "on", x$n - x$df.residual, "Df\n")
    }

    invisible(x)
 }
