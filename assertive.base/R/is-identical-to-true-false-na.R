#' @rdname Truth
#' @export
is_identical_to_false <- function(x, allow_attributes = FALSE, 
  .xname = get_name_in_parent(x))
{
  if(allow_attributes) 
  {
    x <- strip_attributes(x)
  }
  if(!identical(FALSE, x)) 
  {
    msg <- gettextf(
      "%s is not identical to FALSE.", 
      .xname, 
      domain = "R-assertive.base"
    )
    return(false(msg))
  }
  TRUE
}                  

#' @rdname Truth
#' @export
is_identical_to_na <- function(x, allow_attributes = FALSE, 
  .xname = get_name_in_parent(x))
{
  if(allow_attributes) 
  {
    x <- strip_attributes(x)
  }
  if(!identical(NA, x) && 
     !identical(NA_real_, x) && 
     !identical(NA_character_, x) && 
     !identical(NA_integer_, x) && 
     !identical(NA_complex_, x))
  {
    msg <- gettextf(
      "%s is not identical to NA.", 
      .xname, 
      domain = "R-assertive.base"
    )
    return(false(msg))
  }
  TRUE
}

#' @rdname Truth
#' @export
is_identical_to_true <- function(x, allow_attributes = FALSE, 
  .xname = get_name_in_parent(x))
{
  if(allow_attributes) 
  {
    x <- strip_attributes(x)
  }
  if(!identical(TRUE, x))
  {
    msg <- gettextf(
      "%s is not identical to TRUE.", 
      .xname, 
      domain = "R-assertive.base"
    )
    return(false(msg))
  }
  TRUE
}
