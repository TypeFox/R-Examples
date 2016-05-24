#' @title Count and Percentage
#'
#' @description A function for calculating and formatting counts and
#' percentages.
#'
#' @details
#' Default behavior will return the count of successes and the percentage as "N
#' (pp%).  If there are missing values in the input vector the omission of such
#' can be controlled by setting \code{na.rm = TRUE}.  In this case, the number
#' of non-missing values will be reported by default.  Omission of the
#' non-missing values can be controlled by setting \code{show_denom = "never"}.
#'
#' The function n_perc0 uses a set of default arguments which may be
#' advantageous for use in building tables. 
#' 
#'
#' @param x a 0:1 or boolean vector
#' @param digits digits to the right of the decimal point to return in the
#' percentage estimate.
#' @param na_rm if true, omit NA values
#' @param show_denom defaults to "ifNA".  Other options are "always" or "never".
#' @param show_symbol if TRUE (default) the percent symbol is shown, else it is supressed.
#' @param markup latex or markdown
#'
#' @return a character vector of the formatted values
#'
#' @examples
#' 
#' n_perc(c(0, 1,1, 1, 0, 0), show_denom = "always")
#' n_perc(c(0, 1,1, 1, 0, 0, NA), na_rm = TRUE)
#' 
#' n_perc(mtcars$cyl == 6)
#'
#' set.seed(42)
#' x <- rbinom(4269, 1, 0.314)
#' n_perc(x)
#' n_perc(x, show_denom = "always")
#' n_perc(x, show_symbol = FALSE)
#'
#' # n_perc0 examples
#' n_perc0(c(0, 1,1, 1, 0, 0))
#' n_perc0(mtcars$cyl == 6)
#'
#' @rdname n_perc
#' @export   
n_perc <- function(x, 
                   digits = getOption("qwraps2_frmt_digits", 2), 
                   na_rm = FALSE, 
                   show_denom = "ifNA", 
                   show_symbol = TRUE,
                   markup = getOption("qwraps2_markup", "latex")) { 
  d <- sum(!is.na(x))
  n <- sum(x, na.rm = na_rm)
  p <- frmt(100 * n/d, digits)

  if (show_denom == "never") { 
    rtn <- paste0(frmt(as.integer(n)), " (", p, "%)")
  } else { 
    if (show_denom =="always" | any(is.na(x))) { 
      rtn <- paste0(frmt(as.integer(n)), "/", frmt(as.integer(d)), " (", p, "%)")
    } else { 
      rtn <- paste0(frmt(as.integer(n)), " (", p, "%)")
    }
  }

  if (!show_symbol) { 
    rtn <- gsub("%", "", rtn)
  }


  if (markup == "latex") { 
    rtn <- gsub("%", "\\\\%", rtn)
  } 

  return(rtn)
}

#' @rdname n_perc
#' @export
perc_n <- function(x, 
                   digits = getOption("qwraps2_frmt_digits", 2), 
                   na_rm = FALSE, 
                   show_denom = "ifNA", 
                   markup = getOption("qwraps2_markup", "latex")) { 
  d <- sum(!is.na(x))
  n <- sum(x, na.rm = na_rm)
  p <- frmt(100 * n/d, digits)

  if (show_denom == "never") { 
    rtn <- paste0(frmt(as.integer(n)), " (", p, "%)")
  } else { 
    if (show_denom =="always" | any(is.na(x))) { 
      rtn <- paste0(frmt(as.integer(n)), "/", frmt(as.integer(d)), " (", p, "%)")
    } else { 
      rtn <- paste0(frmt(as.integer(n)), " (", p, "%)")
    }
  }

  if (markup == "latex") { 
    rtn <- gsub("%", "\\\\%", rtn)
  } 

  return(rtn)
}

#' @rdname n_perc
#' @export   
n_perc0 <- function(x, 
                   digits = 0,
                   na_rm = FALSE, 
                   show_denom = "never", 
                   show_symbol = FALSE,
                   markup = getOption("qwraps2_markup", "latex")) { 
  n_perc(x, digits, na_rm, show_denom, show_symbol, markup)
}

