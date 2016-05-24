################################################################################
##
##   R package clusrank by Mei-Ling Ting Lee, Jun Yan, and Yujing Jiang
##   Copyright (C) 2015
##
##   This file is part of the R package clusrank.
##
##   The R package clusrank is free software: you can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package clusrank is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package clusrank. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################


#'Print Function for Nonparametric Test for Clustered Data
#'
#'Printing objects of class "\code{ctest}", by simple \link{print}
#'methods.
#'
#'@author Yujing Jiang
#'
#'@param x object of class "\code{ctest}".
#'@param digits number of significant digits to be used.
#'@param prefix string, passed to \link{strwrap}
#'for displaying the method component of the htest object.
#'@param ... further arguments to be passed to or from methods.
#'
#' @return the argument \code{x}, invisibly, as for all \link{print}
#' methods.
#' @examples
#' data(crd)
#'clt <- cluswilcox.test(z ~ group(group) + cluster(id), data = crd)
#'print(clt, digits = 2)
#'@export

print.ctest <- function (x, digits = getOption("digits"), prefix = "\t", ...)
{
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  cat("data:  ", x$data.name, "\n", sep = "")
  out <- character()
  if (!is.null(x$rstatistic)){
     cat(paste(names(x$rstatistic), "=", format(signif(x$rstatistic,
                                                                max(1L, digits - 2L)))))
  cat("\n")
  }
  if (!is.null(x$erstatistic)){
    cat( paste(names(x$erstatistic), "=", format(signif(x$erstatistic,
                                                                max(1L, digits - 2L)))))

  cat("\n")}

  if (!is.null(x$vrstatistic)){
    cat(paste(names(x$vrstatistic), "=", format(signif(x$vrstatistic,
                                                                 max(1L, digits - 2L)))))
  cat("\n")}

  out <- character(0)
  if (!is.null(x$statistic))
    out <- c(out, paste(names(x$statistic), "=", format(signif(x$statistic,
                                                               max(1L, digits - 2L)))))
  if (!is.null(x$parameter))
    out <- c(out, paste(names(x$parameter), "=", format(signif(x$parameter,
                                                               max(1L, digits - 2L)))))
  if (!is.null(x$p.value)) {
    fp <- format.pval(x$p.value, digits = max(1L, digits -
                                                3L))
    out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) ==
                                       "<") fp else paste("=", fp)))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  out <- character()
  if(!is.null(x$n))
    out <- c(out, paste(names(x$n), "=", format(signif(x$n, max(1L, digits - 2L)))))
  if(!is.null(x$cn))
    out <- c(out, paste(names(x$cn), "=", format(signif(x$cn,
                                                       max(1L, digits - 2L)))))
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")

  if(!is.null(x$adjusted) && x$adjusted == TRUE)
    cat("The signed rank test statistics is adjusted since the data is unbalanced.")
  if(!is.null(x$balance) && x$balance == FALSE)
    cat("The data is unbalanced")

    cat("\n")
  if (!is.null(x$alternative)) {
    cat("alternative hypothesis: ")
    if (!is.null(x$null.value)) {
      if (length(x$null.value) == 1L) {
        alt.char <- switch(x$alternative, two.sided = "not equal to",
                           less = "less than", greater = "greater than")
        cat("true ", names(x$null.value), " is ", alt.char,
            " ", x$null.value, "\n", sep = "")
      }
      else {
        cat(x$alternative, "\nnull values:\n", sep = "")
        print(x$null.value, digits = digits, ...)
      }
    }
    else cat(x$alternative, "\n", sep = "")
  }
  if (!is.null(x$conf.int)) {
    cat(format(100 * attr(x$conf.int, "conf.level")), " percent confidence interval:\n",
        " ", paste(format(c(x$conf.int[1L], x$conf.int[2L])),
                   collapse = " "), "\n", sep = "")
  }
  if (!is.null(x$estimate)) {
    cat("sample estimates:\n")
    print(x$estimate, digits = digits, ...)
  }
  cat("\n")
  invisible(x)
}
