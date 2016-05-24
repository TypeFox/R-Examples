#' @rdname is_r
#' @export
is_r_alpha <- function()
{
  if(version$status != "alpha")
  {
    return(not_this_build("alpha"))
  }
  TRUE
}

#' @rdname is_r
#' @export
is_r_beta <- function()
{
  if(version$status != "beta")
  {
    return(not_this_build("beta"))
  }
  TRUE
}

#' @rdname is_r
#' @export
is_r_devel <- function()
{
  if(version$status != "Under development (unstable)")
  {
    return(not_this_build("development"))
  }
  TRUE
}

#' @rdname is_r
#' @export
is_r_patched <- function()
{
  if(version$status != "Patched")
  {
    return(not_this_build("patched"))
  }
  TRUE
}

#' @rdname is_r
#' @export
is_r_release_candidate <- function()
{
  if(version$status != "RC")
  {
    return(not_this_build("release candidate"))
  }
  TRUE
}

#' @rdname is_r
#' @export
is_r_release <- function()
{
  if(nzchar(version$status))
  {
    return(not_this_build("release"))
  }
  TRUE
}

#' @rdname is_r
#' @export
is_r_revised <- function()
{
  if(version$status != "Revised")
  {
    return(not_this_build("revised"))
  }
  TRUE
}

#' @rdname is_r
#' @export
is_r_stable <- function()
{
  .Deprecated("is_r_release")
  is_r_release()
}

#' Failure for bad build
#' 
#' Wrapper to \code{false} for failure messages when the OS is not as 
#' expected.
#' @param status A string giving the name of the build status that was desired.
#' @return A string showing the actual build status.
#' @seealso \code{\link[base]{.Platform}} and \code{\link[base]{Sys.info}}
#' @examples
#' \donttest{
#' assertive.reflection:::not_this_build("stable")
#' assertive.reflection:::not_this_build("development")
#' }
#' @noRd
not_this_build <- function(status)
{
  reported_status <- clean_status_string()
  false(
    gettextf(
      "R's build type is %s, not %s.",
      reported_status,
      status
    )
  ) 
}

clean_status_string <- function(status = version$status)
{
  switch(
    status,
    Patched                        = "patched",
    "Under development (unstable)" = "development",
    alpha                          = "alpha",
    beta                           = "beta",
    RC                             = "release candidate",
    Revised                        = "revised",
    "release"
  )
}
