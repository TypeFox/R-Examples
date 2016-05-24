#' Does the string match a pattern?
#' 
#' Checks to see if in the input matches a regular expression or fixed character
#' pattern.
#' @param x string
#' @param pattern pattern
#' @param opts_fixed Passed to \code{\link[stringi]{stri_detect_fixed}}.
#' @param opts_regex Passed to \code{\link[stringi]{stri_detect_regex}}.
#' @param na_ignore should NAs be ignored or not?
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @seealso \code{\link[stringi]{stri_detect}}, on which these functions are 
#' based.
#' @examples
#' # Is it safe to eat oysters?
#' is_matching_fixed(month.name, "r")
#' 
#' # Sometimes it is easier to specify the negative match.
#' is_matching_regex(LETTERS, "[^AEIOU]")
#' is_not_matching_regex(LETTERS, "[AEIOU]")
#' 
#' # Matching is vectorized over both x and pattern
#' (pi_digits <- strsplit(format(pi, digits = 17), "")[[1]])
#' is_matching_regex(pi_digits, c("[13]", "[59]"))
#' 
#' assert_any_are_matching_regex(pi_digits, c("[13]", "[59]"))
#' 
#' # These checks should fail
#' assertive.base::dont_stop({
#'   assert_all_are_matching_regex(pi_digits, c("[13]", "[59]"))
#' })
#' @author Aditya Bhagwat
#' @rdname is_matching_fixed
#' @importFrom stringi stri_detect_fixed
#' @importFrom assertive.base false get_name_in_parent
#' @export
is_matching_fixed <- function(x, pattern, opts_fixed = NULL, 
                              .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x)
    {
      ok <- stringi::stri_detect_fixed(x, pattern, opts_fixed = opts_fixed)
      set_cause(ok, gettextf("does not match '%s'", pattern))
    },
    x
  )
}

#' @rdname is_matching_fixed
#' @importFrom stringi stri_detect_fixed
#' @importFrom assertive.base false get_name_in_parent
#' @export
is_not_matching_fixed <- function(x, pattern, opts_fixed = NULL, 
                              .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x)
    {
      ok <- !stringi::stri_detect_fixed(x, pattern, opts_fixed = opts_fixed)
      set_cause(ok, gettextf("matches '%s'", pattern))
    },
    x
  )
}


#' @rdname is_matching_fixed
#' @importFrom stringi stri_detect_regex
#' @importFrom assertive.base false get_name_in_parent
#' @export
is_matching_regex <- function(x, pattern, opts_regex = NULL, 
                              .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x)
    {
      ok <- stringi::stri_detect_regex(x, pattern, opts_regex = opts_regex)
      set_cause(ok, gettextf("does not match '%s'", pattern))
    },
    x
  )
}

#' @rdname is_matching_fixed
#' @importFrom stringi stri_detect_regex
#' @importFrom assertive.base false get_name_in_parent
#' @export
is_not_matching_regex <- function(x, pattern, opts_regex = NULL, 
                                  .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x)
    {
      ok <- !stringi::stri_detect_regex(x, pattern, opts_regex = opts_regex)
      set_cause(ok, gettextf("matches '%s'", pattern))
    },
    x
  )
}