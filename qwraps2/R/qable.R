#' @title Qable: an extended verion of knitr::kable
#'
#' @description Create a simple table via kable with row groups and rownames
#' similar to those of \code{hmisc::latex} or \code{htmlTable::htmlTable}.
#'
#' @details
#' TO DO
#'
#' @seealso
#' hmisc::latex, htmlTable::htmlTable
#'
#' @param x \code{matrix} or \code{data.frame} to be turned into a qable
#' @param rgroup a named numeric vector with the name of the row group and the
#' number of rows within the group.  \code{sum(rowgroup) == nrow(x)}.
#' @param rnames a character vector of the row names
#' @param cnames column names
#' @param markup the markup language to use
#' @param ... additional arguments passed to \code{knitr::kable}
#'
#' @return a character vector of the formatted numbers
#'
#' @examples
#' library(dplyr)
#' 
#' this_summary <- function(.data) { 
#'   summarize(.data, 
#'             qwraps2::frmt(min(mpg)), 
#'             qwraps2::frmt(median(mpg)),
#'             qwraps2::frmt(max(mpg)), 
#'             qwraps2::frmt(min(hp)), 
#'             qwraps2::frmt(max(hp)), 
#'             qwraps2::frmt(mean(wt)))
#' }
#' 
#' mtcars$cyl_factor <- factor(mtcars$cyl, levels = c(4, 6, 8))
#' 
#' tab <- cbind(mtcars %>% this_summary %>% t,
#'              mtcars %>% group_by(cyl_factor) %>% this_summary %>% t %>% {.[-1, ]})
#' 
#' rwgrp <- c("Miles Per Gallon" = 3, "Horse Power" = 2, "Weight" = 1)
#' rwnms <- c("Min MPG", "Median MPG", "Max MPG", "Min HP", "Max HP", "Mean Weight")
#' cnms  <- c("All mtcars", paste(levels(mtcars$cyl_factor), "Cyl"))
#' 
#' qable(tab, rwgrp, rwnms, cnms, markup = "latex")
#' qable(tab, rwgrp, rwnms, cnms, markup = "markdown")
#'
#' @export   
#' @rdname qable
qable <- function(x, rgroup, rnames, cnames, markup = getOption("qwraps2_markup", "latex"), ...) { 

  rg_idx <- cumsum(c(1, 1 + rgroup[-length(rgroup)]))

  if (markup == "latex") { 
    xmat <- matrix("~", nrow = nrow(x) + length(rgroup), ncol = 1 + ncol(x))
    xmat[rg_idx, 1] <- paste0("\\bf{", names(rgroup), "}")
    xmat[-rg_idx, 1] <- paste("~~", rnames)
  } else if (markup == "markdown") { 
    xmat <- matrix("&nbsp;", nrow = nrow(x) + length(rgroup), ncol = 1 + ncol(x))
    xmat[rg_idx, 1] <- paste0("**", names(rgroup), "**")
    xmat[-rg_idx, 1] <- paste("&nbsp;&nbsp;", rnames)
  } else {
    stop("markup is either 'latex' or 'markdown'")
  }

  xmat[-rg_idx, -1] <- x

  knitr::kable(xmat, row.names =  FALSE, col.names = c("", cnames), ...)
}

