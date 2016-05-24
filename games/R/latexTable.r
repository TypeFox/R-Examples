##
## INPUT:
## x: character vector
##
## RETURN:
## x, with common special characters in LaTeX converted to their escape
## sequences
## 
latexEsc <- function(x)
{
    x <- gsub("{", "\\{", x, fixed = TRUE)
    x <- gsub("}", "\\}", x, fixed = TRUE)
    x <- gsub("_", "\\_", x, fixed = TRUE)
    x <- gsub("#", "\\#", x, fixed = TRUE)
    x <- gsub("$", "\\$", x, fixed = TRUE)
    x <- gsub("%", "\\%", x, fixed = TRUE)
    x <- gsub("^", "\\^", x, fixed = TRUE)
    x <- gsub("&", "\\&", x, fixed = TRUE)
    x <- gsub("~", "\\textasciitilde{}", x, fixed = TRUE)
    return(x)
}

##' LaTeX table for strategic models
##' 
##' Makes a LaTeX table of strategic model results.
##'
##' \code{latexTable} prints LaTeX code for the presentation of results from a
##' strategic model.  Each row contains one regressor, and each column contains
##' one of the utility (or variance term) equations in the model.  For example,
##' in a model fit with \code{\link{egame12}}, the four columns will be u11,
##' u13, u14, and u24 respectively.  Each cell contains the estimated parameter,
##' atop its standard error in parentheses.  Cells that are empty because a
##' regressor is not included in a particular equation are filled with the
##' string specified in the option \code{blankfill}.  Signorino and Tarar (2006,
##' p. 593) contains a table of this variety.
##'
##' The table generated depends on the \pkg{multirow} package for LaTeX, so
##' make sure to include \code{\\usepackage{multirow}} in the preamble of your
##' document.
##' 
##' The \code{digits} option does not yet work seamlessly; you may have to
##' resort to trial and error.
##' @param x a fitted model of class \code{game}.
##' @param digits number of digits to print.
##' @param scientific logical or integer value to control use of scientific
##' notation.  See \code{\link{format}}.
##' @param blankfill text to fill empty cells (those where a certain variable
##' did not enter the given equation).
##' @param math.style.negative whether negative signs should be "math style" or
##' plain hyphens.  Defaults to \code{TRUE}.
##' @param file file to save the output in.  Defaults to \code{""}, which prints
##' the table to the R console.
##' @param floatplace where to place the table float; e.g., for
##' \code{\\begin\{table\}[htp]}, use \code{floatplace = "htp"}.
##' @param caption caption to use (none if \code{NULL})
##' @param label reference label to use (none if \code{NULL})
##' @param rowsep amount of space (in points) to put between rows.
##' @param useboot logical: use bootstrap estimates (if available) to calculate
##' standard errors?
##' @return \code{x}, invisibly.
##' @export
##' @references Curtis S. Signorino and Ahmer Tarar.  2006.  "A Unified Theory
##' and Test of Extended Immediate Deterrence."  \emph{American Journal of
##' Political Science} 50(3):586--605.
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
##' @example inst/examples/latexTable.r
latexTable <- function(x, digits = max(3, getOption("digits") - 2), scientific =
                       NA, blankfill = "", math.style.negative = TRUE, file =
                       "", floatplace = "htbp", caption = NULL, label = NULL,
                       rowsep = 2, useboot = TRUE)
{
    ## this whole function should be rewritten using writeLines; this is my hack
    ## until then
    lcat <- function(...) cat(..., file = file, sep = "", append = TRUE)

    ## get coefficient names, values, and standard errors
    n <- names(coef(x))
    cf <- summary(x, useboot = useboot)$coefficients[, 1:2]
    cf <- rbind(cf, c(sum(x$log.likelihood), 0))

    ## convert coefficients and standard errors to character strings with the
    ## proper number of digits and the correct negative sign
    cf <- format(cf, digits = digits, trim = TRUE, scientific = scientific)
    if (math.style.negative)
        cf <- gsub("-", "$-$", cf)

    ## retrieve equation names and the variable names associated with each.
    ## this is currently done using some character-string hacks; see
    ## 'print.game' comments (in games.r) for an idea on how to do this better.
    eqNames <- x$equations[attr(x$equations, "hasColon")]
    varNames <- sapply(sapply(strsplit(n, ":"), '[', -1), paste, collapse = ":")
    varNames <- unique(varNames[nchar(varNames) > 0])
    if ("(Intercept)" %in% varNames)  # force intercept to be first
        varNames <- c("(Intercept)", varNames[varNames != "(Intercept)"])
    otherNames <- x$equations[!attr(x$equations, "hasColon") &
                              x$equations %in% n[!x$fixed]]

    ## top matter
    lcat("\n%% latex table generated in R ", as.character(getRversion()),
         " by games package\n")
    lcat("%% ", date(), "\n")
    lcat("%% remember to include \\usepackage{multirow} in your preamble\n\n")
    lcat("\\begin{table}[", floatplace, "]\n")
    lcat("\\begin{center}\n")
    lcat("\\begin{tabular}{",
         paste(c("l", rep("c", length(eqNames))), collapse = ""), "}\n")
    lcat("\\hline\n")
    lcat(paste(c("", latexEsc(eqNames)), collapse = " & "), " \\\\\n")
    lcat("\\hline\n")

    ## each row of the table: a variable name
    for (i in varNames) {
        lcat("\\multirow{2}{*}{", latexEsc(i), "} & ")
        ## each column: an equation
        ## first inner loop: coefficient values
        for (J in 1:length(eqNames)) {
            j <- eqNames[J]
            ji <- paste(j, i, sep = ":")
            if (ji %in% n) {
                lcat(cf[ji, 1])
            } else {
                lcat("\\multirow{2}{*}{", blankfill, "}")
            }
            if (J < length(eqNames))
                lcat(" & ")
        }
        lcat(" \\\\\n & ")
        ## second inner loop: standard errors
        for (J in 1:length(eqNames)) {
            j <- eqNames[J]
            ji <- paste(j, i, sep = ":")
            if (ji %in% n)
                lcat("(", cf[ji, 2], ")")
            if (J < length(eqNames))
                lcat(" & ")
        }
        lcat(" \\\\[", rowsep, "pt]\n")

    }

    ## 'otherNames' are the equations estimated only with an intercept, and the
    ## variance terms in the ultimatum model.  it may be good to make it
    ## optional to treat these equations differently
    if (length(otherNames) > 0) {
        lcat("\\hline\n")
        for (i in otherNames) {
            lcat("\\multirow{2}{*}{", i, "} & ", cf[i, 1], " \\\\\n & (",
                 cf[i, 2], ") \\\\[", rowsep, "pt]\n")
        }
    }

    ## bottom matter
    lcat("\\hline \\hline\n")
    lcat("Log-likelihood & ", cf[nrow(cf), 1], " \\\\\n $N$ & ",
         nrow(x$model), "\\\\\n")
    lcat("\\hline\n")
    
    lcat("\\end{tabular}\n")
    lcat("\\end{center}\n")

    if (!is.null(caption))
        lcat("\\caption{", caption, "}\n")
    if (!is.null(label))
        lcat("\\label{", label, "}\n")
    lcat("\\end{table}\n")

    invisible(x)
}
