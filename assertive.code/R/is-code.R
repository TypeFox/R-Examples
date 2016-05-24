#' Does the current call have an argument?
#'
#' Checks to see if the current call has an argument with 
#' the name given in the input.
#'
#' @param x Argument to check. 
#' @param fn Function to find the argument in.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{has_arg} reimplements \code{\link[methods]{hasArg}}, 
#' letting you choose the function to search in, and providing more
#' information on failure.  
#' @note \code{has_arg} is for interactive use and takes an unquoted argument 
#' name; \code{has_arg_} is for programmatic use and takes a string naming a 
#' argument.
#' @seealso \code{\link[methods]{hasArg}}.
#' @examples
#' has_arg(x, mean.default)
#' has_arg(y, mean.default)   
#' f <- function(...) has_arg(z)   
#' f(z = 123)
#' f(123)
#' @importFrom methods formalArgs
#' @export
has_arg <- function(x, fn = sys.function(sys.parent()))
{
  x <- get_name_in_parent(x)
  has_arg_(x, fn)
}

#' @rdname has_arg
#' @export
has_arg_ <- function(x, fn = sys.function(sys.parent()))
{
  formal_args_of_fn <- formalArgs(fn)
  if(!x %in% formal_args_of_fn)
  {                             
    fn_name <- get_name_in_parent(fn)
    fail <- false(
      gettext("%s is not an argument of %s."), 
      sQuote(x), 
      sQuote(fn_name)
    )
    if("..." %in% formal_args_of_fn)
    {
      dots_call <- eval(quote(substitute(list(...))), sys.parent())
      if(!x %in% names(dots_call))
      {
        return(fail)
      }
    } else
    {
      return(fail)
    }
  }
  TRUE
}

#' Is the binding of a variable locked?
#' 
#' Check to see if the binding of a variable is locked (that is, it has been 
#' made read-only).
#' @param x Input to check. (Unlike \code{bindingIsLocked}, you can pass the
#' variable itself, rather than a string naming that variable.)
#' @param env Environment to check where binding had been locked.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{TRUE} or \code{FALSE}, depending upon whether or not the 
#' binding is locked in the specified environment.
#' \code{assert_is_binding_locked} returns nothing but throws an error if 
#' the corresponding \code{is_*} function returns \code{FALSE}.
#' @note The environment is guessed as follows: The name of \code{x} is 
#' determined via \code{get_name_in_parent}.  Then find is called, 
#' @seealso \code{\link[base]{bindingIsLocked}}, which this wraps, 
#' \code{\link[utils]{find}} for how the environment is guessed.  If this returns
#' a single environment, that is used.  Otherwise the parent environment is 
#' used (as determined with \code{\link[base]{parent.frame}}).
#' @examples
#' is_binding_locked(a_non_existent_variable)
#' e <- new.env()
#' e$x <- 1:10
#' is_binding_locked(x, e)
#' lockBinding("x", e)
#' is_binding_locked(x, e)
#' unlockBinding("x", e)
#' is_binding_locked(x, e)
#' @importFrom utils find
#' @importFrom assertive.properties is_scalar
#' @importFrom assertive.types assert_is_environment
#' @export
is_binding_locked <- function(x, env = if(is_scalar(e <- find(.xname))) as.environment(e) else parent.frame(), .xname = get_name_in_parent(x))
{
  assert_is_environment(env)
  .xname <- force(.xname)
  env <- force(env)
  if(!exists(.xname, env, inherits = FALSE))
  {
    return(
      false(
        gettext("%s does not exist in %s."), 
        .xname, 
        format(env)
      )
    )
  }
  if(!bindingIsLocked(.xname, env))
  {
    return(
      false(
        gettext("%s is not locked (read-only) in %s."), 
        .xname, 
        format(env)
      )
    )
  }
  TRUE
}

#' Is the input function being debugged?
#'
#' Checks to see if the input DLL (a.k.a. shared object) is loaded.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_debugged} wraps \code{\link[base]{isdebugged}}, providing 
#' more information on failure.  \code{assert_is_debugged} returns nothing but
#' throws an error if \code{is_debugged} returns \code{FALSE}.
#' @seealso \code{\link[base]{isdebugged}}.
#' @importFrom assertive.base use_first
#' @export
is_debugged <- function(x, .xname = get_name_in_parent(x))
{
  # isdebugged accepts x as either a function or a string
  if(!is.function(x))
  {
    x <- coerce_to(use_first(x), "character", .xname)
  }
  if(!isdebugged(x))
  {
    return(
      false(
        gettext("%s is not being debugged."), 
        .xname
      )
    )
  }
  TRUE
}

#' Does the code run without throwing an error?
#' 
#' Call the code inside a try block and report if an error was thrown.
#' 
#' @param x Code to check.
#' @note Note that this has the side effect of running the code contained in
#' \code{x}.
#' @return \code{TRUE} if the code runs without throwing an error.  The result
#' of running the code is contained in an attribute named \code{"result"}.
#' @export
is_error_free <- function(x)
{
  res <- try(x, silent = TRUE)
  if(inherits(res, "try-error"))
  {
    return(false(attr(res, "condition")$message))
  }
  ok <- TRUE
  attr(ok, "result") <- res
  ok
}

#' Does the variable exist?
#'
#' Checks to see if the input variables exist.
#'
#' @param x Input to check.
#' @param envir Passed to \code{exists}.
#' @param inherits Passed to \code{exists}.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_existing} is a vectorized wrapper to \code{exists}, 
#' providing more information on failure (and with a simplified interface).  
#' The \code{assert_*} functions return nothing but throw an error if 
#' \code{is_existing} returns \code{FALSE}.
#' @seealso \code{\link[base]{exists}}.
#' @examples
#' e <- new.env()
#' e$x <- 1
#' e$y <- 2
#' assert_all_are_existing(c("x", "y"), envir = e)
#' #These examples should fail.
#' assertive.base::dont_stop(assert_all_are_existing(c("x", "z"), envir = e))
#' @importFrom assertive.properties is_empty
#' @importFrom assertive.base bapply
#' @export
is_existing <- function(
  x, 
  envir = parent.frame(), 
  inherits = TRUE, 
  .xname = get_name_in_parent(x)
)
{
  x <- coerce_to(x, "character", .xname)
  if(is_empty(x)) return(logical(0))
  if(length(x) > 1L)
  {
    return(bapply(
      x, 
      is_existing,
      envir    = envir,
      inherits = inherits
    ))
  }
  if(!exists(
    x, 
    envir    = envir,
    inherits = inherits
  ))
  {
    return(false(gettext("%s does not exist."), .xname))
  }
  TRUE
}

#' Is suitable to be used as an if condition
#' 
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_if_condition} returns \code{TRUE} if the input is 
#' scalar \code{TRUE} or \code{FALSE}.
#' @note \code{if} will try to do the right thing if you pass it a number
#' or a string, but this function assumes you want to do the right thing
#' and pass either \code{TRUE} or \code{FALSE}, maybe with some attributes.
#' @examples
#' is_if_condition(TRUE)
#' is_if_condition(FALSE)
#' is_if_condition(NA)
#' is_if_condition(c(TRUE, FALSE))
#' is_if_condition("the truth")
#' # You can pass a number as a logical condition, but you shouldn't,
#' # so the next line returns FALSE.
#' is_if_condition(1)
#' assertive.base::dont_stop(assert_is_if_condition(raw(1)))
#' @importFrom assertive.properties is_scalar
#' @importFrom assertive.types is_logical
#' @importFrom assertive.base is_na
#' @export
is_if_condition <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_logical(x, .xname)))
  {
    return(ok)
  }
  if(!(ok <- is_scalar(x, "length", .xname)))
  {
    return(ok)
  }
  if(is_na(x))
  {
    return(false("%s is NA.", .xname))
  }
  TRUE
}

#' Is the input a symbol in a loaded DLL?
#'
#' Checks to see if the input DLL (a.k.a. shared object) is loaded.
#'
#' @param x A string naming the symbol in a DLL to check.
#' @param PACKAGE A string naming an R package to restrict the search to, or 
#' \code{""} to check all packages. Passed to \code{is.loaded}.
#' @param type A string naming the type of external code call to restrict the
#' search to, or \code{""} to check all type. Passed to \code{is.loaded}.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_loaded} wraps \code{\link[base]{is.loaded}}, providing more 
#' information on failure.  \code{assert_is_loaded} returns nothing but
#' throws an error if \code{is_loaded} returns \code{FALSE}.
#' @seealso \code{\link[base]{is.loaded}}.
#' @export
is_loaded <- function(x, PACKAGE = "", type = c("", "C", "Fortran", "Call", "External"), 
  .xname = get_name_in_parent(x))
{
  type <- match.arg(type)
  if(nzchar(PACKAGE))
  {
    if(is.null(getLoadedDLLs()[[PACKAGE]]))
    {
      return(false(gettext("The DLL %s is not loaded."), PACKAGE))
    }
    routines <- getDLLRegisteredRoutines(PACKAGE)
    type <- if(type == "") 
    {
      c(".C", ".Fortran", ".Call", ".External")
    } else
    {
      paste0(".", type)
    }
    routine_names <- unlist(lapply(type, function(x) names(routines[[x]])))
    if(!x %in% routine_names)
    {
      return(
        false(
          gettext("The routine %s is not registered with the DLL %s."), 
          .xname, 
          PACKAGE
        )
      )
    }
  }
  if(!is.loaded(x, PACKAGE = PACKAGE, type = type))
  {
    return(false(gettext("The symbol %s is not loaded."), .xname))
  }
  TRUE
}

#' Is the input valid R code?
#'
#' Check to see if the input is a valid (parseable) R code.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{TRUE} if the input string is valid R code.
#' @examples
#' is_valid_r_code("x <- 1 + sqrt(pi)")
#' is_valid_r_code("x <- ")
#' is_valid_r_code("<- 1 + sqrt(pi)")
#' @seealso \code{\link[base]{parse}}
#' @importFrom assertive.base use_first
#' @export
is_valid_r_code <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  x <- use_first(x)
  ok <- is_error_free(parse(text = x))
  if(!ok)
  {
    return(false(
      gettext("%s is not valid R code. %s."), 
      .xname, 
      cause(ok)
    ))
  }
  TRUE
}

#' Is the string a valid variable name?
#'
#' Checks strings to see if they are valid variable names.
#'
#' @param x Input to check.
#' @param allow_reserved If \code{TRUE} then "..." and "..1", "..2", etc. 
#' are considered valid.
#' @param allow_duplicates Deprecated and ignored.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return The \code{assert_*} functions return nothing but throw an error if the 
#' corresponding \code{is_*} function returns \code{FALSE}.
#' @seealso \code{\link{make.names}}.
#' @examples
#' make_random_string <- function(n)
#' {
#'   paste0(sample(letters, n, replace = TRUE), collapse = "")
#' }
#' long <- c(make_random_string(10000), make_random_string(10001))
#' x <- c("x", "y_y0.Y", ".", "x y", "...", "..1", long)
#' unname(is_valid_variable_name(x))
#' unname(is_valid_variable_name(x, allow_reserved = FALSE))
#' #These examples should fail.
#' assertive.base::dont_stop(
#'   assert_all_are_valid_variable_names(c("...", "..1"), allow_reserved = FALSE)
#' )
#' @references
#' \url{http://4dpiecharts.com/2011/07/04/testing-for-valid-variable-names/}
#' @importFrom assertive.base set_cause
#' @export
is_valid_variable_name <- function(x, allow_reserved = TRUE, 
  allow_duplicates)
{
  if(!missing(allow_duplicates))
  {
    .Deprecated(
      msg = "The 'allow_duplicates' argument is deprecated and will be ignored."
    )
  }
  x <- coerce_to(x, "character", get_name_in_parent(x))
  
  #is name too long?
  max_name_length <- if(getRversion() < "2.13.0") 256L else 10000L
  ok <- short_enough <- nchar(x) <= max_name_length
  
  not_missing_and_ok <- !is.na(ok) & ok
  
  #is it a reserved variable, i.e.
  #an ellipsis or two dots then a number?
  not_reserved <- rep.int(TRUE, length(x))
  if(!allow_reserved)
  {
    rx <- "^\\.{2}[[:digit:]]+$"
    ok[not_missing_and_ok] <- not_reserved[not_missing_and_ok] <- 
      x[not_missing_and_ok] != "..." & !grepl(rx, x[not_missing_and_ok])
  } 
  
  #are names valid (and maybe unique)
  not_missing_and_ok <- !is.na(ok) & ok
  
  ok[not_missing_and_ok] <- x[not_missing_and_ok] == 
    make.names(x[not_missing_and_ok])
  
  names(ok) <- x
  set_cause(
    ok, 
    ifelse(
      short_enough, 
      ifelse(not_reserved, "bad format", "reserved"),
      "too long"  
    )
  )
}

