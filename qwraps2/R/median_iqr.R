#' @title Median and Inner Quartile Range
#'
#' @description A function for calculating and formatting the median and 
#' inner quartile range of a data vector.
#' #' @details
#' Given a numeric vector, \code{median_iqr} will return a character string with
#' the median and IQR.  Formating of the output will be extended in
#' future versions.
#'
#' @param x a numeric vector
#' @param digits digits to the right of the decimal point to return.
#' @param na_rm if true, omit NA values
#' @param show_n defaults to "ifNA".  Other options are "always" or "never".
#' @param markup latex or markdown
#'
#' @return a character vector of the formatted values
#'
#' @examples
#' set.seed(42)
#' x <- rnorm(1000, 3, 4)
#' median(x)
#' quantile(x, probs = c(1, 3)/4)
#' median_iqr(x)
#' median_iqr(x, show_n = "always")
#' 
#' x[187] <- NA
#' # median_iqr(x) ## Will error
#' median_iqr(x, na_rm = TRUE)
#'
#'
#' @export   
median_iqr <- function(x, 
                       digits = getOption("qwraps2_frmt_digits", 2), 
                       na_rm = FALSE, 
                       show_n = "ifNA", 
                       markup = getOption("qwraps2_markup", "latex")) { 
  n <- sum(!is.na(x))
  m <- stats::median(x, na.rm = na_rm)
  qs <- stats::quantile(x, probs = c(1, 3) / 4, na.rm = na_rm)

  rtn <- paste0(qwraps2::frmt(m, digits), " (",
                qwraps2::frmt(qs[1L], digits), ", ", 
                qwraps2::frmt(qs[2L], digits), ")")

  if (all(!(show_n %in% c("ifNA", "always", "never")))) { 
    warning("'show_n' should be in c('ifNA', 'always', 'never').  Setting to 'ifNA'.")
    show_n <- "ifNA"
  }

  if (show_n == "always" | (show_n == "ifNA" & any(is.na(x)))) { 
    rtn <- paste0(qwraps2::frmt(as.integer(n), digits), "; ", rtn)
  } 

  return(rtn)
}
