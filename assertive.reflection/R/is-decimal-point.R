#' What does the current locale specify for the decimal point?
#' 
#' Does the current locale specify a comma or a period for the decimal point?
#' 
#' @param dp Character to be used as a decimal point.
#' @param type Decimal point for numbers or money?
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_comma_for_decimal_point} returns \code{TRUE} when the 
#' current locale uses a comma for a decimal place, as determined by 
#' \code{Sys.localeconv}.  Similarly, \code{is_period_for_decimal_point} returns 
#' \code{TRUE} when the current locale uses a period (a.k.a. full stop) for a 
#' decimal place.  If R has been compiled without support for locales, then the 
#' value will always be \code{NA}.
#' @references \url{http://www.cplusplus.com/reference/clocale/lconv/}
#' @seealso \code{\link[base]{Sys.localeconv}}
#' @examples
#' # Current settings:
#' is_comma_for_decimal_point()
#' is_comma_for_decimal_point("money")
#' # Or equivalently:
#' is_period_for_decimal_point()
#' is_period_for_decimal_point("money")
#' # A useful guess for reading in files:
#' read_csv <- if(is_comma_for_decimal_point()) read.csv else read.csv2 
#' \dontrun{
#' # Force locale and test (may require admin rights)
#' current_locale <- sys_get_locale()
#' a_period_locale <- if(is_windows()) 
#' {
#'   "English_United Kingdom.1252"
#' } else if(is_mac()) 
#' {
#'   "en_GB"
#' } else if(is_linux()) 
#' {
#'   "en_GB.utf8"
#' } else 
#' {
#'   "en"
#' }
#' sys_set_locale(LC_ALL = a_period_locale)
#' assert_is_period_for_decimal_point()
#' a_comma_locale <- if(is_windows())
#' {
#'   "French_France.1252"
#' } else if(is_mac()) 
#' {
#'   "fr_FR"
#' } else if(is_linux()) 
#' {
#'   "fr_FR.utf8" 
#' } else 
#' {
#'   "fr"
#' }
#' sys_set_locale(LC_ALL = a_comma_locale)
#' assert_is_comma_for_decimal_point()
#' suppressWarnings(sys_set_locale(l = current_locale))
#' }
#' @importFrom assertive.base na
is_xxx_for_decimal_point <- function(dp, type = c("numbers", "money"))
{
  locale_conventions <- Sys.localeconv()
  if(is.null(locale_conventions))
  {
    return(na(gettext("R has been compiled without support for locales.")))
  }
  type <- match.arg(type)
  element <- switch(type, numbers = "decimal_point", money = "mon_decimal_point")
  if(locale_conventions[element] != dp)
  {
    if(!nzchar(locale_conventions[element]))
    {
      return(
        na(
          gettextf(
            "The locale convention for a (%s) decimal point has not been defined.",
            switch(type, numbers = "numeric", money = "monetary")
          )
        )
      )
    }
    return(
      false(
        gettextf(
          "The locale convention is to use a '%s' for a (%s) decimal point.", 
          locale_conventions[element],
          switch(type, numbers = "numeric", money = "monetary")
        )          
      )
    )
  }
  TRUE
}

#' @rdname is_xxx_for_decimal_point
#' @export
is_comma_for_decimal_point <- function(type = c("numbers", "money"))
{
  is_xxx_for_decimal_point(",", type)
}

#' @rdname is_xxx_for_decimal_point
#' @export
is_period_for_decimal_point <- function(type = c("numbers", "money"))
{
  is_xxx_for_decimal_point(".", type)
}
