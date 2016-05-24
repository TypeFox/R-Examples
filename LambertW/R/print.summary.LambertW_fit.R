#' @rdname LambertW_fit-methods
#' @description
#' \code{print.summary.LambertW_fit} tries to be smart about formatting the
#' coefficients, standard errors, etc. and also displays "significance stars" 
#' (like in the output of \code{summary.lm}).
#' @export
print.summary.LambertW_fit <- function(x, ...) {
  cat("Call: ")
  print(x$call)
  cat("Estimation method: ")
  cat(x$method)
  cat("\n")
  cat("Input distribution: ")
  cat(x$distname)
  cat("\n")
  
  cat("\n Parameter estimates:\n")
  if (x$method == "IGMM") {
    cat(" Note: standard errors are only asymptotic, simulation based.\n")
    cat(" If you want more accurate estimates see ?bootstrap .")
  }
  
  printCoefmat(x$coefmat, signif.stars = TRUE)
  if (x$type == "s") {
    cat("-------------------------------------------------------------- \n")
    if (!any(x$distname == c("exp", "chi", "gamma", "F"))) {
      M <- rbind(x$support, x$data.range)
      colnames(M) <- c("a", "b")
      rownames(M) <- c("Support", "Data range")
      print(M)
      
      cat("\n p_m1 = P(non-principal branch affects solution): ")
      cat(x$p_m1)
      cat("\n")
    }
  } else if (x$type == "hh") {
    cat("-------------------------------------------------------------- \n")
    
    cat("\np-value for 'H_0:symmetric' versus 'H_1:skewed': ")
    cat(x$symmetry.p.value)
    cat("\n")
  }
  
  if (x$type != "hh" && (x$distname == "normal" || x$method == "IGMM")) {
    cat("-------------------------------------------------------------- \n")
    if (x$method == "IGMM") {
      x$theta <- tau2theta(x$tau, beta = x$tau[c("mu_x", "sigma_x")])
    }
    cat("\nGiven these input parameter estimates the moments of the output random variable are \n",
        " (assuming Gaussian input): \n ")
    moments.y <- round(unlist(mLambertW(theta = x$theta, distname = "normal")), 2)
    cat("mu_y = ", moments.y[1], 
        "; sigma_y = ", moments.y[2], 
        "; skewness = ", moments.y[3], 
        "; kurtosis = ", moments.y[4], ".\n", sep = "")
    cat("\n")
  }
} 