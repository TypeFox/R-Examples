#' @rdname is_empty_character
#' @importFrom assertive.types is_a_string
#' @export
is_an_empty_string <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_a_string(x))) return(ok)
  if(!is_empty_character(x)) 
  {
    return(
      false(
        gettext("%s contains characters."), 
        .xname
      )
    )
  }
  TRUE
}

#' @rdname is_empty_character
#' @importFrom assertive.types is_a_string
#' @export
is_a_non_empty_string <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_a_string(x))) return(ok)
  if(!is_non_empty_character(x))
  {
    return(
      false(
        gettext("%s has no characters."), 
        .xname
      )
    )
  }
  TRUE
}

#' @rdname is_empty_character
#' @importFrom assertive.types is_a_string
#' @export
is_a_missing_or_empty_string <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_a_string(x, .xname))) return(ok)
  if(!is_missing_or_empty_character(x)) 
  {
    return(
      false(
        gettext("%s is not missing and contains characters."),
        .xname
      )
    )
  }
  TRUE
}

#' @rdname is_empty_character
#' @importFrom assertive.types is_a_string
#' @export
is_a_non_missing_nor_empty_string <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_a_string(x))) return(ok)
  if(!is_non_missing_nor_empty_character(x)) 
  {
    return(
      false(
        gettext("%s is missing or has no characters."),
        .xname
      )
    )
  }
  TRUE
}
