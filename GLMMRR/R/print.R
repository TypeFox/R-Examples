#' Print RRglmGOF values
#'
#' @param x
#' an object of class RRglmGOF.
#' @param digits
#' minimal number of \emph{significant digits} (default: 3).
#' @param ...
#' further arguments passed to or from other methods.
#' @method print RRglmGOF
#' @export
print.RRglmGOF <- function(x, digits = 3, ...)
{
  cat("\n")
  cat("GLMMRR - Binary Randomized Response Data")
  cat("\n\n")
  cat("Goodness-of-Fit Testing\n")
  cat("Response variable:\t\t", x$vars[1], "\n")
  cat("Predictor(s):\t\t\t", x$vars[2:length(x$vars)], "\n")
  cat("Entries dataset:\t\t", x$n, "\n")

  if(any(c(x$pearson$do, x$deviance$do, x$hlemeshow$do)))
  {
    # Create empty data frame
    df.output <- data.frame("Statistic" = numeric(0), "P.value" = numeric(0), "df" = numeric(0), "Groups" = numeric(0))

    if(x$pearson$do)
    {
      df.pearsonOutput <- data.frame("Statistic" = x$pearson$stat, "P.value" = x$pearson$pvalue, "df" = x$pearson$df, "Groups" = x$pearson$nGroup)
      rownames(df.pearsonOutput) <- "Pearson"
      df.output <- rbind(df.output, df.pearsonOutput)
    }

    if(x$deviance$do)
    {
      df.devianceOutput <- data.frame("Statistic" = x$deviance$stat, "P.value" = x$deviance$pvalue, "df" = x$deviance$df, "Groups" = x$deviance$nGroup)
      rownames(df.devianceOutput) <- "Deviance"
      df.output <- rbind(df.output, df.devianceOutput)
    }

    if(x$hlemeshow$do)
    {
      df.hlemeshowOutput <- data.frame("Statistic" = x$hlemeshow$stat, "P.value" = x$hlemeshow$pvalue, "df" = x$hlemeshow$df, "Groups" = x$hlemeshow$nGroup)
      rownames(df.hlemeshowOutput) <- "Hosmer-Lemeshow"
      df.output <- rbind(df.output, df.hlemeshowOutput)
    }


    cat("---------------------------------------------------------\n")
    cat("Summary: \n\n")
    print(df.output, digits = digits, right = FALSE)
    cat("\n")
    if(x$hlemeshow$do)
    {
      cat("---------------------------------------------------------\n")
      cat("Hosmer-Lemeshow: \n\n")
      print(x$hlemeshow$overview, digits = digits, right = FALSE)
      cat("\n")
    }


  }
}
