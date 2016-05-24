#' Does the input contain empty or missing strings?
#' 
#' Checks for empty or missing strings.
#' @param x A character vector.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return The \code{is_*} functions return logical vectors for strings which
#' are (non) empty or missing, and the \code{assert_*} functions throw errors
#' on failure.
#' @note In R, \code{NA_character_} is considered to be a non-empty string
#' (at least by \code{\link[base]{nzchar}}), which is why many functions are
#' needed to to clarify the situation.
#' @seealso \code{\link[assertive.types]{is_character}}, 
#' \code{\link[base]{nzchar}}
#' @examples 
#' # These functions return a vector:
#' x <- c("", "a", NA)
#' is_empty_character(x)
#' is_non_empty_character(x)
#' is_missing_or_empty_character(x)
#' is_non_missing_nor_empty_character(x)
#' 
#' # These functions return a single value:
#' is_an_empty_string("")
#' is_an_empty_string("a")
#' is_an_empty_string(NA_character_)
#' 
#' is_a_non_empty_string("")
#' is_a_non_empty_string("a")
#' is_a_non_empty_string(NA_character_)
#' 
#' is_a_missing_or_empty_string("")
#' is_a_missing_or_empty_string("a")
#' is_a_missing_or_empty_string(NA_character_)
#' 
#' is_a_non_missing_nor_empty_string("")
#' is_a_non_missing_nor_empty_string("a")
#' is_a_non_missing_nor_empty_string(NA_character_)
#' @importFrom assertive.base coerce_to
#' @importFrom assertive.base call_and_name
#' @importFrom assertive.base set_cause
#' @export
is_empty_character <- function(x, .xname = get_name_in_parent(x))
{ 
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x)
    {
      ok <- !nzchar(x)
      set_cause(ok, ifelse(is.na(x), "missing", "nonempty"))
    },
    x
  )
}

#' @rdname is_empty_character
#' @importFrom assertive.base coerce_to
#' @importFrom assertive.base call_and_name
#' @importFrom assertive.base set_cause
#' @export
is_non_empty_character <- function(x, .xname = get_name_in_parent(x))
{ 
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x)
    {
      ok <- nzchar(x)
      set_cause(ok, "empty")
    },
    x
  )
}

#' @rdname is_empty_character
#' @importFrom assertive.base coerce_to
#' @importFrom assertive.base is_na
#' @importFrom assertive.base set_cause
#' @export
is_missing_or_empty_character <- function(x, .xname = get_name_in_parent(x))
{ 
  x <- coerce_to(x, "character", .xname)
  ok <- !nzchar(x) | is_na(x)
  set_cause(ok, "nonempty")
}

#' @rdname is_empty_character
#' @importFrom assertive.base coerce_to
#' @importFrom assertive.base is_na
#' @importFrom assertive.base set_cause
#' @export
is_non_missing_nor_empty_character <- function(x, .xname = get_name_in_parent(x))
{ 
  x <- coerce_to(x, "character", .xname)
  ok <- nzchar(x) & !is_na(x)
  set_cause(ok, ifelse(is.na(x), "missing", "empty"))
}

#' @rdname is_empty_character
#' @export
is_not_missing_nor_empty_character <- function(x)
{
  .Deprecated("is_non_missing_nor_empty_character")
  is_non_missing_nor_empty_character(x)
}


