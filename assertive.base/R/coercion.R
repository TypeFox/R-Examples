#' Alternative version of is
#' 
#' If a function named \code{is.class} exists, call \code{is.class(x)}.
#' If not, call \code{is(x, class)}.
#' @param x Input to check.
#' @param class Target class that \code{x} maybe belong to.
#' @param .xname Not intended to be used directly.
#' @return \code{TRUE} if x belongs to the class and \code{FALSE} 
#' otherwise. 
#' @seealso \code{\link[methods]{is}}, and 
#' \code{\link[assertive.types]{assert_is_all_of}} for the corresponding assert fns.
#' @examples
#' is2(1:5, "character")
#' is2(matrix(1:5), "character")
#' is2(1:5, c("character", "list", "numeric"))
#' @importFrom methods is
#' @export
is2 <- function(x, class, .xname = get_name_in_parent(x))
{    
  # Can't use is_empty in next line because that function calls this one.
  if(length(class) == 0L) stop("You must provide a class.")
  if(length(class) > 1L) 
  {
    return(
      set_cause(
        bapply(class, function(cl) is2(x, cl, "")),
        rep.int(type_description(x), length(class))
      )
    )
  }
  ok <- tryCatch(
    {
      is.class <- match.fun(paste0("is.", class))
      is.class(x)
    },
    error = function(e)
    {
      is(x, class) 
    }
  )
  if(!ok)
  {
    return(
      false(
        "%s is not of type '%s'; it has %s.", 
        .xname, 
        class, 
        type_description(x)
      )
    )
  }
  TRUE
}

#' Coerce variable to a different class
#'
#' Coerce the input to a different class, with a warning.  More reliable then 
#' \code{\link[methods]{as}}, and supports coercion to multiple classes.
#'
#' @param x Input to coerce.
#' @param target_class The desired class of x.  Multiple values allowed (see 
#' note).
#' @param .xname Not intended to be used directly.
#' @return The input \code{x} after attempted coercion to the target class.
#' @note If x does not already have the target class, a warning is given
#' before coercion.  
#' The function will try and convert the \code{x} to each of the classes given
#' in \code{target_class}, in order, until it succeeds or runs out of classes
#' to try.  It will first try and convert \code{x} using a dedicated 
#' \code{as.target_class} function if that exists.  If it does not exist, or 
#' throws an error then \code{coerce_to} will try to use 
#' \code{as(x, target_class)}.
#' @seealso \code{\link[methods]{is}} and \code{\link[methods]{as}}.
#' @examples
#' # Numbers can be coerced to characters but not to calls.
#' dont_stop(coerce_to(1:5, c("call", "character")))
#' @importFrom methods as
#' @export
coerce_to <- function(x, target_class, .xname = get_name_in_parent(x))
{
  # Can't use is_empty in next line because that function calls this one.
  if(length(target_class) == 0L) 
  {
    stop("You must provide a class.")
  }
  if(!is.character(target_class))
  {
    stop("target_class should be a character vector.")
  }
  for(this_class in target_class)
  {
    if(!is2(x, this_class))
    {
      warning(
        sprintf(
          "Coercing %s to class %s.", 
          .xname,
          sQuote(this_class)
        ),
        call. = FALSE
      )
    }
    tryCatch(
      {
        as.this_class <- match.fun(paste0("as.", this_class))
        return(as.this_class(x))
      },
      error = function(e)
      {
        # as.this_class doesn't exist; try as(, "this_class") instead
        tryCatch(
          return(as(x, this_class)),
          error = function(e)
          {
            # Can't coerce to this class; warn and move to next class
            warning(
              sprintf(
                "%s cannot be coerced to type %s.", 
                .xname,
                sQuote(this_class)
              ), 
              call. = FALSE
            )
          }
        )
      }
    )
  }
  # Nothing worked; throw an error
  stop(
    sprintf(
      "%s cannot be coerced to any of these types: %s.", 
      .xname,
      toString(sQuote(target_class))
    )
  )
}

#' Describe the type of object
#' 
#' Get the class or mode (for arrays).
#' @param x A variable.
#' @return A string.
type_description <- function(x)
{
  if(is.array(x))
  {
    sprintf("mode '%s'", mode(x))
  } else
  {
    sprintf("class '%s'", toString(class(x)))
  }
}
