#' Get or set the system locale
#'
#' Wrappers to \code{Sys.getlocale} and \code{Sys.setlocale} for getting and
#' setting the system locale.
#'
#' @param simplify If \code{TRUE}, the locale settings are returned as a vector,
#' otherwise, a list.
#' @param ... Name-value pairs of locale categories to set.
#' @param l A list, as an alternative method of passing local categories to set.
#' @return A named list or vector giving the system locale names. 
#' \code{sys_set_locale} invisibly returns the locale settings *before* making 
#' changes (like \code{setwd} and \code{options} do).
#' @examples
#' (current_locale <- sys_get_locale())
#' \dontrun{
#' english <- if(is_windows()) "English" 
#'   else if(is_mac()) "en_GB" 
#'   else if(is_linux()) "en_GB.utf8" 
#'   else "en"
#' sys_set_locale(LC_MONETARY = english)
#' sys_get_locale()
#' sys_set_locale(l = current_locale)  #restore everything
#' }
#' @seealso \code{\link[base]{Sys.getlocale}}.
#' @export
sys_get_locale <- function(simplify = FALSE)
{
  locale <- Sys.getlocale()
  if(locale == "C")
  {
    categories <- locale_categories(FALSE)
    values <- lapply(categories, function(x) "C")
  } else
  {
    splitter <- if(is_windows() || is_linux()) ";" else "/"
    locale <- strsplit(locale, splitter)[[1]]
    locale <- strsplit(locale, "=")
    categories <- vapply(
      locale,
      function(x) x[1],
      character(1)
    )
    values <- lapply(
      locale,
      function(x) x[2]
    )
  }
  
  names(values) <- categories
  if(simplify) unlist(values) else values
}

#' @rdname sys_get_locale
#' @importFrom assertive.base merge_dots_with_list
#' @export
sys_set_locale <- function(..., l = list())
{
  old_locale <- sys_get_locale()
  values <- merge_dots_with_list(..., l = l)
  categories <- names(values)
  categories <- match.arg(
    categories,
    locale_categories(),
    several.ok = TRUE
  )
  
  for(i in seq_along(values))
  {
    Sys.setlocale(categories[i], values[[i]])
  }
  invisible(old_locale)
}

#' Allowed locale categories.
#'
#' The categories of locale that can be gotten/set.
#'
#' @param include_all If \code{TRUE}, the value \code{LC_ALL} is included.
#' @param include_unix If \code{TRUE}, the extra unix-only values are included.
#' @return A character vector of locale categories.
#' @seealso \code{\link{sys_get_locale}}.
locale_categories <- function(include_all = TRUE, include_unix = is_unix())
{
  allowed_categories <- c(
    if(include_all) "ALL",
    "COLLATE", "CTYPE", "MONETARY", "NUMERIC", "TIME",
    if(include_unix) c("MESSAGES", "PAPER", "MEASUREMENT")
  )
  paste0("LC_", allowed_categories)
}
