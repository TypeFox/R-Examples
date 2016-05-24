#' @title Format Wrappers
#'
#' @description Functions for formating Numeric values for consistent display
#' in reports.
#'
#' @details
#' `frmt` is really just a wrapper for the \code{formatC}.
#'
#' `frmtp` formats P-values per journal
#' requirements.  As I work on papers aimed at different journals, the
#' formatting functions will be extended to match.
#'
#' Default settings are controlled through the function arguments but should be
#' set via \code{options()}.
#' 
#' Default settings report the P-value exactly if P >
#' \code{getOptions("qwraps2_frmtp_digits", 4)} and reports 
#' P < \code{10^-(getOptions("qwraps2_frmtp_digits", 2))} otherwise.  By the
#' leading zero is controlled via 
#' \code{getOptions("qwraps2_frmtp_leading0", TRUE)} 
#' and a upper or lower case P is controlled by
#' \code{getOptions("qwraps2_frmtp_case", "upper")}.  These options are ignored 
#' if \code{style != "default"}.
#'
#' Journals with predefined P-value formatting are noted in the
#' \pkg{\link{qwraps2}} documentation.
#'
#' `frmtci` takes a \code{matrix}, or \code{data.frame}, with a point estimate
#' and the lcl and ucl and formats a string for reporting.  est (lcl, ucl) is
#' the default.  The confidence level can be added to the string, e.g., "est
#' (95% CI: lcl, ucl)" with \code{show_level = TRUE} or by setting your own
#' format.
#'
#' `frmtcip` expects four values, est, lcl, ucl, and p-value.  The resulting
#' sting will be of the form "est (lcl, ucl; p-value)".
#' 
#'
#' @seealso
#' \code{\link{formatC}}
#'
#' @param x a vector of numbers or a numeric matrix to format.
#' @param digits number of digits, including trailing zeros, to the right of the
#' decimal point.  This option is ignored if \code{is.integer(x) == TRUE)}.
#'
#' @return a character vector of the formatted numbers
#'
#' @examples
#' 
#' # Formatting numbers
#' integers <- c(1234L, 9861230L)
#' numbers  <- c(1234,  9861230)
#' frmt(integers)  # no decimal point
#' frmt(numbers)   # decimal point and zeros to the right
#' 
#' numbers <- c(0.1234, 0.1, 1234.4321, 0.365, 0.375)
#' frmt(numbers)
#'
#' # Formatting p-values
#' ps <- c(0.2, 0.001, 0.00092, 0.047, 0.034781, 0.0000872, 0.787, 0.05, 0.043)
#' # LaTeX is the default markup language
#' cbind("raw"      = ps, 
#'       "default"  = frmtp(ps), 
#'       "3lower"   = frmtp(ps, digits = 3, case = "lower"),
#'       "PediDent" = frmtp(ps, style = "pediatric_dentistry"))
#'
#' # Using markdown
#' cbind("raw"      = ps, 
#'       "default"  = frmtp(ps, markup = "markdown"), 
#'       "3lower"   = frmtp(ps, digits = 3, case = "lower", markup = "markdown"),
#'       "PediDent" = frmtp(ps, style = "pediatric_dentistry", markup = "markdown"))
#'
#' @export   
#' @rdname frmt
frmt <- function(x, digits = getOption("qwraps2_frmt_digits", 2)) { 
  sapply(x, 
         function(xx) { 
           if (is.integer(xx)) { 
             formatC(xx, format = "d", big.mark = ",")
           } else { 
             formatC(xx, digits = digits, format = "f", big.mark = ",")
           }
         })
}


#' @export   
#' @rdname frmt
#' @param style a character string indicating a specific journal requirements
#' for p-value formatting.  
#' @param markup a character string indicating if the output should be latex or
#' markup.
#' @param case a character string indicating if the output should be upper case
#' or lower case.
#' @param leading0 boolean, whether or not the p-value should be reported as
#' 0.0123 (TRUE, default), or .0123 (FALSE).
frmtp <- function(x,
                  style    = getOption("qwraps2_journal", "default"), 
                  digits   = getOption("qwraps2_frmtp_digits", 4), 
                  markup   = getOption("qwraps2_markup", "latex"),
                  case     = getOption("qwraps2_frmtp_case", "upper"), 
                  leading0 = getOption("qwraps2_frmtp_leading0", TRUE)) {  

  rtn <- 
    switch(style,
           default               = frmtp_default(x, digits, case, leading0), 
           obstetrics_gynecology = frmtp_obstetrics_gynecology(x), 
           pediatric_dentistry   = frmtp_pediatric_dentistry(x))

  if (markup == "latex") { 
    rtn <- paste0("$", rtn, "$")
  } else if (markup == "markdown") {
    rtn <- gsub("(P|p)", "\\*\\1\\*", rtn)
  }
  return(rtn)
}

frmtp_default <- function(x, digits, case, leading0) { 
  p_cutoff <- 10^-digits
  
  sapply(x, 
         function(xx) { 
           if (xx < p_cutoff) { 
             p_val <- paste("P <", formatC(p_cutoff, format = "g"))
           } else {
             p_val <- paste0("P = ", formatC(xx, digits = digits, format = "f")) 
           }

           if (case == "lower") { 
             p_val <- tolower(p_val)
           }

           if (!leading0) { 
             p_val <- gsub("0\\.", "\\.", p_val)
           }

           p_val
         })
}

frmtp_pediatric_dentistry <- function(x) { 
  sapply(x, function(xx) {

         if (xx < 0.001) { 
           p_val <- "< .001"
         } else if (xx < 0.01) { 
           p_val <- paste0("= ", formatC(xx, digits = 3, format = "f"))
         } else if (xx < 0.05 & round(xx, 2) == 0.05) { 
           p_val <- paste0("= ", formatC(xx, digits = 3, format = "f"))
         } else { 
           p_val <- paste0("= ", formatC(xx, digits = 2, format = "f"))
                             }
         paste0("P ", gsub("0\\.", "\\.", p_val))
      })
}

frmtp_obstetrics_gynecology <- function(x) { 
  sapply(x, function(xx) {

         if (xx < 0.001) { 
           p_val <- "< .001"
         } else { 
           p_val <- paste0("= ", formatC(xx, digits = 3, format = "f"))
                             }
         paste0("P ", gsub("0\\.", "\\.", p_val))
      })
}

#' @export
#' @rdname frmt
#' @param est the numeric index of the vector element or the matrix column
#' containing the point estimate.
#' @param lcl the numeric index of the vector element or the matrix column
#' containing the lower confidence limit.
#' @param ucl the numeric index of the vector element or the matrix column
#' containing the upper confidence limit.
#' @param format a string with "est" "lcl", and "ucl" to denote the location of
#' the estimate, lower confidence limit, and upper confidence limit for the
#' formated string.  Defaults to "est (lcl, ucl)".
#' @param show_level defualts to FALSE.  If TRUE and \code{format} is the
#' default, then
#' "100*(1-options()$qwraps2_alpha)% CI:" will be placed between the left
#' parenthesis and the lcl.  If set to a string, then the given string will be
#' placed between the left parenthesis and the lcl.  If the \code{format} is not
#' the default, then this argument is ignored.
#' @param ... args passed to frmt
#' @examples
#' # Formatting the point estimate and confidence interval
#' # for a set of three values
#' temp <- c(a = 1.23, b = .32, CC = 1.78)
#' frmtci(temp) 
#' frmtci(temp, show_level = TRUE) 
#' 
#' # note that the show_level will be ignored in the following
#' frmtci(temp, format = "est ***lcl, ucl***", show_level = TRUE)
#' 
#' # show_level as a character
#' frmtci(temp, show_level = "confidence between: ")
#' 
#' # For a matrix: the numbers in this example don't mean anything, but the
#' # formatting should.
#' temp2 <- matrix(rnorm(12), nrow = 4, 
#'                 dimnames = list(c("A", "B", "C", "D"), c("EST", "LOW", "HIGH")))
#' temp2
#' frmtci(temp2)
frmtci <- function(x, est = 1, lcl = 2, ucl = 3, format = "est (lcl, ucl)", show_level = FALSE, ...) {
  UseMethod("frmtci") 
}

#' @export
frmtci.default <- function(x, est = 1, lcl = 2, ucl = 3, format = "est (lcl, ucl)", show_level = FALSE, ...) {

  .est <- qwraps2::frmt(x[est], ...)
  .lcl <- qwraps2::frmt(x[lcl], ...)
  .ucl <- qwraps2::frmt(x[ucl], ...)

  if (format == "est (lcl, ucl)") { 
    out <- sub("ucl", .ucl, sub("lcl", .lcl, sub("est", .est, format))) 
    if (is.character(show_level)) {
      lvl <- paste0("(", show_level)
      out <- sub("\\(", lvl, out) 
    } else if(is.logical(show_level)) { 
      if (show_level) { 
      lvl <- paste0("(", 100 * (1 - getOption("qwraps2_alpha", 0.05)), "% CI: ")
      out <- sub("\\(", lvl, out)
      }
    }
  } else { 
    out <- sub("ucl", .ucl, sub("lcl", .lcl, sub("est", .est, format)))
  } 
  out
}

#' @export
frmtci.matrix <- function(x, est = 1, lcl = 2, ucl = 3, format = "est (lcl, ucl)", show_level = FALSE, ...) {
  apply(x, 1, frmtci.default, est = 1, lcl = lcl, ucl = ucl, format = format, show_level = show_level, ...)
}

#' @export
frmtci.qwraps2_mean_ci <- function(x, est = 1, lcl = 2, ucl = 3, format = "est (lcl, ucl)", show_level = FALSE, ...) {

  .est <- qwraps2::frmt(x[est], ...)
  .lcl <- qwraps2::frmt(x[lcl], ...)
  .ucl <- qwraps2::frmt(x[ucl], ...)

  if (format == "est (lcl, ucl)") { 
    out <- sub("ucl", .ucl, sub("lcl", .lcl, sub("est", .est, format))) 

    if (show_level) { 
      lvl <- paste0("(", 100 * (1 - attr(x, "alpha")), "% CI: ")
      out <- sub("\\(", lvl, out)
    }
  } else { 
    out <- sub("ucl", .ucl, sub("lcl", .lcl, sub("est", .est, format)))
  } 
  out
}
