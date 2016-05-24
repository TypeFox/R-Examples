#' Wrapper to vapply that returns booleans
#' 
#' Wrapper to \code{\link{vapply}} for functions that return a boolean (logical 
#' scalar) value.
#' 
#' @param x A vector (atomic or list).
#' @param predicate A predicate (function that returns a bool) to apply.
#' elementwise to \code{x}.
#' @param ... Passed to \code{vapply}.
#' @return A logical vector.
#' @note \code{USE.NAMES} is set to \code{TRUE}
#' @seealso \code{\link{vapply}}.
#' @export
bapply <- function(x, predicate, ...)
{
  vapply(x, predicate, logical(1L), ..., USE.NAMES = TRUE)
}

#' Call a function, and give the result names.
#'
#' Calls a function, and names the result with the first argument.
#'
#' @param fn A function to call.  See note below.
#' @param x The first input to \code{fn}.
#' @param ... Optional additional inputs to \code{fn}.
#' @return The result of \code{fn(x, ...)}, with names given by the
#' argument \code{x}.
#' @note The function, \code{fn}, should return an object with the 
#' same length as the input \code{x}.  For speed and simplicity, this
#' isn't checked; it is up to the developer of the assertion to make
#' sure that this condition holds.
#' @examples
#' call_and_name(is.finite, c(1, Inf, NA))
#' @seealso \code{\link{cause}} and \code{\link{na}}.
#' @export
call_and_name <- function(fn, x, ...)
{
  y <- fn(x, ...)
  dim(y) <- dim(x)
  names(y) <- to_names(x)
  y
}

# Lots of issues about how best to generate names!
# Original behaviour was to use as.character everywhere, but high precision
# wanted for numbers, so behaviour was changed to use format.  This breaks 
# lots of tests.  Here are the requirements.
# - For atomic vectors, NA values should be given a missing name, not "NA"
#   -> Can't use deparse, format
# - Numbers should be given to high precision (inc complex)
#   -> Can't use as.character on numbers
# - Character vectors shouldn't be quoted
#   -> Can't use deparse
# - Recursive variables should just be a deparse, but exact details not too fussy (too rare)
to_names <- function(x)
{
  if(is.double(x))
  {
    ifelse(is.na(x), NA_real_, sprintf("%.17g", x))
  } else if(is.complex(x))
  {
    ifelse(is.na(x), NA_complex_, sprintf("%.17g+%.17gi", Re(x), Im(x)))
  } else
  {
    as.character(x)
  }
}

#' Run code without stopping
#' 
#' Runs code without stopping for warnings or errors.
#' @param expr Code to execute.
#' @return A list containing the results of evaluating each call in \code{expr}.
#' @note This function is dangerous, since it overrides warnings and errors.
#' Its intended use is for documenting examples of warnings and errors.
#' @seealso \code{\link[base]{warning}} and \code{\link[base]{stop}} for 
#' generating warnings and errors respectively; \code{\link[base]{try}} and
#' \code{\link[base]{conditions}} for handling them.
#' @examples
#' dont_stop({
#'   warning("a warning")
#'   x <- 1
#'   stop("an error")
#'   y <- sqrt(exp(x + 1))
#'   assert_is_identical_to_true(y)
#'   y > 0
#' })
#' @export
dont_stop <- function(expr)
{
  this_env <- sys.frame(sys.nframe())
  
  # Split the expression up into a list of calls
  subbed_expr <- substitute(expr, this_env)
  # Temporarily wrap expr in braces, if it isn't already
  brace <- quote(`{`)
  if(!identical(subbed_expr[[1L]], brace))
  {
    subbed_expr <- c(brace, subbed_expr)
  }
  call_list <- as.list(subbed_expr)[-1L] # -1 to ignore brace again
  names(call_list) <- vapply(
    call_list, 
    function(x) 
    {
      paste0(deparse(x), collapse = "")
    }, 
    character(1)
  )
  
  handler <- function(e)
  {
    e["call"] <- list(NULL) # Override the condition's call
    e
  }
  
  # Evaluate each one in turn
  lapply(
    call_list,
    function(calli)
    {
      tryCatch(eval(calli, this_env), warning = handler, error = handler)
    }
  )
}

#' Get the name of a variable in the parent frame
#'
#' Gets the name of the input in the parent frame.
#'
#' @param x Variable to get the name of.
#' @param escape_percent Logical. If \code{TRUE}, percent signs are doubled, 
#' making the value suitable for use with \code{sprintf} (and hence by 
#' \code{false} and \code{na}).
#' @return A string giving the name of the input in the parent frame.
#' @examples 
#' outside <- 1
#' f <- function(inside, escape_percent) 
#' {
#'   get_name_in_parent(inside, escape_percent)
#' }
#' f(outside, TRUE) 
#' f('10%', TRUE) 
#' f('10%', FALSE)
#' @export
get_name_in_parent <- function(x, escape_percent = TRUE)
{  
  xname <- paste0(
    deparse(
      do.call(
        substitute, 
        list(substitute(x), parent.frame())
      )
    ),
    collapse = ""
  )
  if(escape_percent)
  {
    xname <- gsub("%", "%%", xname)
  }
  xname
}

#' Merge two lists
#' 
#' Merges two lists, taking duplicated elements from the first list.
#' @param x A list.
#' @param y A list.
#' @param warn_on_dupes \code{TRUE} or \code{FALSE}.  Should a warning be given 
#' if both \code{x} and \code{y} have elements with the same name.  See note.
#' @param allow_unnamed_elements \code{TRUE} or \code{FALSE}. Should unnamed
#' elements be allowed?
#' @param ... Ignored.
#' @return A list, combining elements from \code{x} and \code{y}.
#' @note In the event of elements that are duplicated between \code{x} and 
#' \code{y}, the versions from \code{x} are used.
#' @seealso \code{\link{merge_dots_with_list}}, \code{\link[base]{merge}}
#' @examples
#' merge(
#'   list(foo = 1, bar = 2, baz = 3), 
#'   list(foo = 4, baz = 5, quux = 6)
#' )
#' 
#' # If unnamed elements are allowed, they are included at the end
#' merge(
#'   list("a", foo = 1, "b", bar = 2, baz = 3, "c"), 
#'   list(foo = 4, "a", baz = 5, "b", quux = 6, "d"),
#'   allow_unnamed_elements = TRUE
#' )
#' @method merge list
#' @export
merge.list <- function(x, y, warn_on_dupes = TRUE, allow_unnamed_elements = FALSE, ...)
{
  if(length(y) == 0) return(x)
  y <- coerce_to(y, "list", get_name_in_parent(y))
  
  # Get elements without names
  x_is_unnamed <- names_never_null(x) == ""
  y_is_unnamed <- names_never_null(y) == ""
  if(allow_unnamed_elements)
  {
    unnamed_values <- c(x[x_is_unnamed], y[y_is_unnamed])
    x <- x[!x_is_unnamed]
    y <- y[!y_is_unnamed]
  } else # !allow_unnamed_elements
  {
    if(any(x_is_unnamed) || any(y_is_unnamed))
    {
      stop("There are unnamed elements in x or y, but allow_unnamed_elements = FALSE.")
    }
  }
  
  # Now deal with named elements
  all_names <- c(names(x), names(y))
  all_values <- c(x, y)
  if(anyDuplicated(all_names) > 0)
  {
    if(warn_on_dupes)
    {
      warning(
        "Duplicated arguments: ", 
        toString(all_names[duplicated(all_names)])
      )
    }
    all_values <- all_values[!duplicated(all_names)]
  }
  if(allow_unnamed_elements)
  {
    all_values <- c(all_values, unnamed_values)
  }
  all_values
}

#' Merge ellipsis args with a list.
#'
#' Merges variable length ellipsis arguments to a function with a list argument.
#'
#' @param ... Some inputs.
#' @param l A list.
#' @param warn_on_dupes \code{TRUE} or \code{FALSE}.  Should a warning be given 
#' if both \code{x} and \code{y} have elements with the same name.  See note.
#' @param allow_unnamed_elements \code{TRUE} or \code{FALSE}. Should unnamed
#' elements be allowed?
#' @note If any arguments are present in both the \code{...} and \code{l} 
#' arguments, the \code{...} version takes preference, and a warning is thrown.
#' @return A list containing the merged inputs.
#' @seealso \code{\link{merge.list}}, \code{\link[base]{merge}}
#' @examples
#' merge_dots_with_list(
#'   foo = 1, 
#'   bar = 2, 
#'   baz = 3, 
#'   l = list(foo = 4, baz = 5, quux = 6)
#' )
#' @export
merge_dots_with_list <- function(..., l = list(), warn_on_dupes = TRUE, allow_unnamed_elements = FALSE)
{
  dots <- list(...)
  l <- coerce_to(l, "list", get_name_in_parent(l))
  merge(dots, l, warn_on_dupes = warn_on_dupes, allow_unnamed_elements = allow_unnamed_elements)
}

#' Wrap a string in brackets
#'
#' Parenthesise a character vector by wrapping elements in brackets, 
#' dashes or commas.
#' @param x Character vector to wrap in parenthenses.
#' @param type String naming the type of parenthesis.
#' @return A character vector of the input wrapped in parentheses.
#' @note English grammar terminology is awfully confusing.  The verb 'to 
#' parenthesise' means to wrap a phrase in brackets or dashes or commas,
#' thus denoting it as supplementary material that could be left out.
#' A 'parenthesis' as a noun is often used as a synonym for a round bracket.
#' @seealso \code{\link[base]{sQuote}}
#' @examples
#' paste("There were three", parenthesise(3), "mice in the experiment.")
#' paste(
#'   "I love parmos", 
#'   parenthesise("Teesside's finest culinary invention", "en_dashes"), 
#'   "but they are sure to give me heart disease."
#' )
#' parenthesise(letters[1:5], "curly")
#' paste0(
#'   "The R language", 
#'   parenthesise("an offshoot of S and Scheme", "commas"), 
#'   "is quite good for data analysis."
#' )
#' @export
parenthesize <- function(x, 
  type = c("round_brackets", "square_brackets", "curly_brackets", "angle_brackets", "chevrons", "hyphens", "en_dashes", "em_dashes", "commas")) 
{
  type <- match.arg(type)
  x <- coerce_to(x, "character", get_name_in_parent(x))
  before <- switch(
    type,
    round_brackets  = "(",
    square_brackets = "[",
    curly_brackets  = "{",
    angle_brackets  = "<",
    chevrons        = "\u3008",
    hyphens         = "- ",
    en_dashes       = "\u2013 ",
    em_dashes       = "\u2014",
    commas          = ", "
  )
  after <- switch(
    type,
    round_brackets  = ")",
    square_brackets = "]",
    curly_brackets  = "}",
    angle_brackets  = ">",
    chevrons        = "\u3009",
    hyphens         = " -",
    en_dashes       = " \u2013",
    em_dashes       = "\u2014",
    commas          = ", "
  )
  paste0(before, x, after)
}

#' @rdname parenthesize
#' @export
parenthesise <- parenthesize

#' Strip all attributes from a variable
#'
#' Strips all the attributes from a variable.
#'
#' @param x Input to strip.
#' @return \code{x}, without attributes.
#' @examples
#' x <- structure(c(foo = 1, bar = 2), some_attr = 3)
#' x2 <- strip_attributes(x)
#' attributes(x)
#' attributes(x2)
#' @export
strip_attributes <- function(x)
{
  attributes(x) <- NULL
  x
}
 
#' Only use the first element of a vector
#'
#' If the input is not scalar, then only the first element is returned, 
#' with a warning.
#'
#' @param x Input that should be scalar.
#' @param indexer Either double indexing, \code{"[["} (the default) or
#' single indexing \code{"["}.
#' @param .xname Not intended to be used directly.
#' @return If \code{x} is scalar, it is returned unchanged, otherwise
#' only the first element is returned, with a warning.
#' @examples 
#' dont_stop(use_first(1:5))
#' @export
use_first <- function(x, indexer = c("[[", "["), .xname = get_name_in_parent(x))
{
  len_x <- length(x)
  # Can't use assert_is_non_empty, is_scalar in next lines because those 
  # functions calls this one.
  if(len_x == 0L)
  {
    stop(sprintf("%s has length 0.", .xname))
  }
  if(len_x == 1L)
  {
    return(x)
  }
  indexer <- match.fun(match.arg(indexer))
  warning(
    sprintf("Only the first value of %s will be used.", .xname),
    call. = FALSE
  )
  indexer(x, 1L)
}
