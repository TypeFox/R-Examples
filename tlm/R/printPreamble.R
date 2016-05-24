printPreamble <-
function(x)
 {
  family <- family(x$model)$family
  modelname <- switch(family,
                      "gaussian" = "Linear",
                      "binomial" = "Logistic",
                      "poisson" = "Poisson")
  modelname <- paste(modelname, "regression") 
  ytrans <- switch(family,
                   "gaussian" = "none",
                   "binomial" = "logit link for logistic regression",
                   "poisson" = "log link for Poisson regression")
  if (family == "gaussian" && x$ypow != 1)
   {
    if (x$ypow == 0) ytrans <- "log" else ytrans <- paste("power, exponent =", attr(x, "ypowlabel"))
   }
  if (x$xpow == 1)
   xtrans <- "none" else
    {
     if (x$xpow == 0) xtrans <- "log" else xtrans <- paste("power, exponent =", attr(x, "xpowlabel"))
    }
  cat("\n")
  if (family == "gaussian" && x$ypow == 1 & x$xpow == 1)
   {
   	cat(modelname, "fitted model \n")
    cat(rep("-", nchar(modelname) + 13, sep = ""), sep = "")
    cat("\n")
    } else {
    cat(modelname, "fitted model in the transformed space\n")
    cat(rep("-", nchar(modelname) + 38, sep = ""), sep = "")
    cat("\n\n")
    cat("Transformations:\n")
    if ((family == "gaussian" && x$ypow != 1) || family != "gaussian")
     cat("   In the response variable:", ytrans, "\n")
    if (x$xpow != 1)
    cat("   In the explanatory variable:", xtrans, "\n")	
   }
 }
