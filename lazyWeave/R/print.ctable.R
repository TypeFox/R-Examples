#' @name WritePrintCtable
#' @export 
#' @method print ctable
#' 
#' @title Write and Print Comparison Tables
#' @description Print Comparisons to the console or write the table to a file.
#' 
#' @param x A ctable object to be printed or written
#' @param max.rows The maximum number of rows to be displayed on a page.  
#'   Variable groups are forced to be kept together.
#' @param round The number of decimal places to be displayed for numeric values.
#' @param percent Toggles if percentages or proportions are printed for 
#'   categorical values
#' @param quartile Toggles if quartiles or min and max are printed for numeric
#'   values associated with the \code{median} argument in \code{conttable}
#' @param cwidth A vector giving the width of each column in the body of the 
#'   table.  The number of columns is not always obvious, and this vector is
#'   not recycled.  If the length is inappropriate, a warning message will
#'   be printed indicating the correct number of columns.
#' @param caption The name of the table.
#' @param footnote A footnote for the table.
#' @param byVarN Toggles if the N per group in \code{byVar} are printed in the
#'   column headings.
#' @param size A character string denoting the size of the text for the 
#'   table.  This must be latex code, for example "\\normalsize" or "\\small."
#'   Remember to use to backslashes!
#' @param descripCombine Toggles if descriptive statistics are combined.  
#'   It is strongly recommended that this be left \code{TRUE} as the appearance
#'   is much better.
#' @param oddsCombine Toggles if odds ratios are combined with the lower up
#'   upper confidence limits.
#' @param markSignificant Toggles if significant results are printed in bold
#'   text.
#' @param statHeader Character string giving the column heading for statistical 
#'   summaries.
#' @param name Toggles if the variable name is printed in the table.
#' @param var.label Toggles if the variable label is printed in the table.
#' @param level Toggles if the variable levels are printed in the table.  This
#'   column is usually needed for categorical variables, and never needed for
#'   numeric variables.
#' @param total Toggles if the totals column is printed in the table
#' @param descriptive Toggles if descriptive statistics are printed.
#' @param missing Toggles if the number of missing values are printed in the 
#'   table.
#' @param missing.perc Toggles if the percentage of missing values is printed 
#'   in the table.
#' @param testStat Toggles if the test statistics are printed.
#' @param odds Toggles if odds ratios are printed.  Only relevant to numeric
#'   variables.
#' @param pval Toggles if the pvalue column is printed.
#' @param oneLine When true, binary variables are printed with only one 
#'   line per variable.  This does not affect the printing of
#'   numeric variable or variables with more than two levels.
#' @param keepVarTogether Determines if variables are kept together when 
#'   splitting a table across mulitple pages.  When a single variable has more 
#'   levels than can remain on one printed page, this must be set to FALSE.
#' @param ... Other arguments to be passed to \code{print} (for 
#'   \code{print.ctable} or \code{lazy.table} (for \code{write.ctable}).  
#'   Currently none are implemented for \code{print}.  In \code{split_ctable}, 
#'   any options in \code{write.ctable} may also be included.
#' @param pvalFormat Character string passed to \code{pvalString} and determines
#'   the pvalue style to be printed.
#' @param pvalArgs A list of additional arguments to be passed to \code{pvalString}
#' @param cat Logical. Determines if the output is returned as a character string
#'   or returned via the \code{cat} function (printed to console).  The default
#'   value is set by \code{options()$lazyWeave_cat}.  This argument allows for
#'   selective override of the default.
#'   
#' @author Benjamin Nutter


print.ctable <- function(x, ...){

  nlev <- nlevels(attributes(x)$byVar)
  lev <- levels(attributes(x)$byVar)

  x <- as.data.frame(x)

  summary.names <- function(r){
    k <- x$type[r]
    if (k == "Bootstrap Mean"){
      n <- c("boot",  "lowerb", "upperb")
      note <- "*a*"
    }
    else if (k == "Parametric Mean"){
      n <- c("mean",   "sd",     "prop")
      note <- "*b*"
    }
    else if (k == "Median"){
      n <- c("median", "p25",    "p75")
      note <- "*c*"
    }
    else{
      n <- c("prop",   "mean",   "sd")
      note <- "*d*"
    }

    n <- as.vector(t(sapply(n, grep, names(x))))
    summ <- x[r, n]
    names(summ) <- paste("stat", rep(LETTERS[1:3], nlev), sep="")
    names(summ) <- paste(names(summ), ".", rep(lev, each=3), sep="")
    if (k %in% c("Bootstrap Mean", "Parametric Mean", "Median") ||
        all(is.na(summ)))
      rownames(summ) <- paste(rownames(summ), note)
    summ
  }
  
  descrip <- round(do.call("rbind", lapply(1:nrow(x), summary.names)), 2)
  out <- cbind(x[, "total", drop=FALSE], descrip,
               x[, c("test.stat", "pvalue")])
  out$test.stat <- round(out$test.stat, 2)
  out$pvalue <- round(out$pvalue, 3)
  out[is.na(out)] <- ""
  rownames(out) <- rownames(descrip)
  print(out)
  cat("\n*a* Mean and Bootstrap Confidence Limits",
      "\n*b* Mean and Standard Deviation",
      "\n*c* Median and Quartiles",
      "\n*d* Percentage\n")
}

