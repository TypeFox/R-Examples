plot.proLik <- function(x, CL = TRUE, alpha = 0.05, type = "l", ...)
{
  if(CL == TRUE & is.null(alpha)) stop("'alpha' must be non-null to include CL in the plot")
  
  plot(x$lambdas ~ x$var.estimates,
     main = deparse(substitute(x)), 
     xlab = x$component, ylab = "LRT statistic", 
     type = type, ...)
     if(CL){  
        chi <- (0.5 * qchisq(alpha, df = 1, lower.tail = FALSE))
        abline(h = chi, lty = "dotted", col = "red", lwd = 2)
    }  
}

