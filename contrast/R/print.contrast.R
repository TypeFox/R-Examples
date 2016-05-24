# This method is used for printing the objects returned by the contrast methods.
# It was copied from the rms package, written by Frank Harrell.
print.contrast <- function(x, X=FALSE, fun=function(u) u, ...)
{
   testLabels <- switch(
      x$model,
      lm =, glm =, lme =, gls = c("t", "Pr(>|t|)"),
      geese = c("Z", "Pr(>|Z|)"))
      
   w <- x[c("Contrast", "SE", "Lower", "Upper", "testStat", "df", "Pvalue")]
   w$testStat <- round(w$testStat, 2)
   w$Pvalue <- round(w$Pvalue, 4)
   no <- names(w)
   no[no == 'SE'] <- 'S.E.'
   no[no == 'testStat'] <- testLabels[1]
   no[no == 'Pvalue'] <- testLabels[2]
   names(w) <- no
   
   cat(x$model, "model parameter contrast\n\n")

   cnames <- x$cnames
   if (length(cnames) == 0)
   {
      cnames <- if (x$nvary) rep('', length(x[[1]])) else as.character(1:length(x[[1]]))
   }
   attr(w, 'row.names') <- cnames
   attr(w, 'class') <- 'data.frame'
   w$Contrast <- fun(w$Contrast)
   w$SE <- fun(w$SE)
   if(x$model != "geese") w$df <- x$df
   w$Lower <- fun(w$Lower)
   w$Upper <- fun(w$Upper)
   print(as.matrix(w), quote=FALSE)
   if (X)
   {
      attr(x$X, "contrasts") <- NULL
      attr(x$X, "assign") <- NULL
      cat('\nContrast coefficients:\n')
      if (is.matrix(x$X)) dimnames(x$X) <- list(cnames, dimnames(x$X)[[2]])
      print(x$X)
   }

   if(x$model == "lm")
   {
      if(x$covType != "const") cat("\nThe", x$covType, "covariance estimator was used.\n")
    }

     
   invisible()
}

