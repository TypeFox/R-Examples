#' Are you running R?
#'
#' Checks to see what type of R you are running.
#'
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_r} wraps \code{is.R}, providing more 
#' information on failure.  \code{is_r_stable}, \code{is_r_patched},
#' \code{is_r_devel}, etc., tell you what type of R build you are 
#' running.  \code{is_architect}, \code{is_rstudio} and \code{is_revo_r} tell
#' you if you are running Architect/StatET, RStudio, or Revolution Analytics'
#' Revolution R build.  \code{is_slave_r} tells you if you are running a slave
#' instance of R (e.g. when building a package with \code{devtools} or using a
#' cluster).
#' The \code{assert_*} functions return nothing but throw an error if 
#' the corresponding \code{is_*} function returns \code{FALSE}.
#' @seealso \code{\link[base]{is.R}}, \code{\link[base]{version}}.
#' @references \url{http://www.revolutionanalytics.com/revolution-r-open}
#' @examples
#' # If this is FALSE, you really need to ditch that old copy of S-PLUS
#' is_r()
#' assertive.base::dont_stop(assert_is_r())
#' # Release, patched, devel, etc.
#' is_r_release()
#' is_r_patched()
#' is_r_devel()
#' is_r_alpha()
#' is_r_beta()
#' is_r_release_candidate()
#' is_r_revised()
#' switch(
#'   version$status,
#'   Patched                        = assert_is_r_patched(),
#'   "Under development (unstable)" = assert_is_r_devel(),
#'   alpha                          = assert_is_r_alpha(),
#'   beta                           = assert_is_r_beta(),
#'   RC                             = assert_is_r_release_candidate(),
#'   Revised                        = assert_is_r_revised(),
#'   assert_is_r_release()
#' )
#' # IDE
#' is_architect()
#' is_revo_r()
#' is_rstudio()
#' @export
is_r <- function()
{
  if(!exists("is.R") || !is.function(is.R) || !is.R())
  {
    return(false("You are not running R."))
  } 
  TRUE
}



#' How is R running?
#' 
#' Tests to see if R is running in batch mode/interactively.
#' 
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_batch_mode} returns \code{TRUE} if R is running in batch 
#' mode.
#' \code{is_interactive} returns \code{TRUE} if R is running interactively.
#' @seealso \code{\link[base]{EnvVar}} and \code{\link[base]{interactive}}.
#' @examples
#' is_batch_mode()
#' is_interactive()
#' is_r_slave()
#' @export
is_batch_mode <- function()
{
  if(is.na(Sys.getenv("R_BATCH", NA)))
  {
    return(false("R is not running in batch mode."))
  }
  TRUE
}

#' @rdname is_batch_mode
#' @export
is_interactive <- function()
{
  if(!interactive())
  {
    return(false("R is not running interactively."))
  }
  TRUE
}

#' @rdname is_batch_mode
#' @export
is_r_slave <- function()
{
  cargs <- commandArgs()
  if(!"--slave" %in% cargs && !all(c("--quiet", "--no-save") %in% cargs))
  {
    return(false("You are not running a slave instance of R."))
  }
  TRUE
}

#' @rdname is_batch_mode
#' @export
is_slave_r <- function()
{
  .Deprecated("is_r_slave")
  is_r_slave()  
}
