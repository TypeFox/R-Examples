#' Can R find tools?
#' 
#' Checks to see if R can see command line tools.
#' @param tools A character vector of tools to look for.
#' @param java_type A string denoting the type of Java to look for (either 
#' 32 or 64 bit).
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return The \code{is_*} functions return \code{TRUE} if the input is 
#' within an interval.  The \code{assert_*} functions return nothing but
#' throw an error if the corresponding \code{is_*} function returns 
#' \code{FALSE}. 
#' @note \code{r_can_compile_code} is a convenience function looking for
#' \code{gcc} and \code{make}.
#' @note \code{r_can_build_translations} is a convenience function looking for
#' \code{gettext} and \code{msgfmt}.
#' @seealso \code{\link[base]{Sys.which}}
#' @examples
#' r_can_find_tools(c("latex", "pdflatex"))
#' r_can_compile_code()
#' r_can_build_translations()
#' r_can_find_java()
#' assertive.base::dont_stop({
#'   assert_r_can_find_tools(c("latex", "pdflatex"))
#'   assert_r_can_compile_code()
#'   assert_r_can_build_translations()
#'   assert_r_can_find_java("64bit")
#' })
#' @export
r_can_find_tools <- function(tools)
{
  paths <- Sys.which(tools)
  not_found <- !nzchar(paths)
  if(any(not_found))
  {
    return(
      false(
        ngettext(
          sum(not_found), 
          "R cannot find the %s tool.", 
          "R cannot find the %s tools."
        ),
        toString(tools[not_found])
      )
    )
  }
  TRUE
}

#' @rdname r_can_find_tools
#' @export
r_can_compile_code <- function()
{
  r_can_find_tools(c("gcc", "make"))
}

#' @rdname r_can_find_tools
#' @export
r_can_build_translations <- function()
{
  r_can_find_tools(c("gettext", "msgfmt"))
}

#' @rdname r_can_find_tools
#' @export
r_can_find_java <- function(java_type = c("same_as_r", "any", "64bit", "32bit"))
{
  java_type <- match.arg(java_type)
  if(!(ok <- r_can_find_tools("java")))
  {
    return(ok) 
  }
  if(java_type == "any")
  {
    return(TRUE)
  }
  if(java_type == "same_as_r")
  {
    java_type <- if(is_64_bit())
    {
      "64bit"
    } else
    {
      "32bit"
    }
  }
  bit_spec <- switch(
    java_type,
    "64bit" = "-d64",
    "32bit" = "-d32"
  )
  cmd <- paste("java -version", bit_spec)
  res <- suppressWarnings(system(cmd, intern = TRUE))
  status <- attr(res, "status")
  if(!is.null(status) && status != 0)
  {
    return(
      false(
        "The %s bit version of Java is not available.", 
        substring(java_type, 1, 2)
      )
    )
  }
  TRUE
}

#' Does R have a capability?
#' 
#' Check to see if R has a specific capability.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return The \code{is_*} functions return \code{TRUE} if R has the capability 
#' and \code{FALSE} (with a cause) otherwise.
#' The \code{assert_*} functions return nothing but throw an error if the 
#' corresponding \code{is_*} function returns \code{FALSE}.
#' @seealso \code{\link[base]{capabilities}}
#' @aliases r_has_capability
#' @examples
#' \dontrun{
#' if(r_has_png_capability())
#' {
#'   png("test.png")
#'   with(cars, plot(speed, dist))
#'   dev.off()
#' } else 
#' {
#'   pdf("test.pdf")
#'   with(cars, plot(speed, dist))
#'   dev.off()
#' }
#' }
#' @export
r_has_jpeg_capability <- function()
{
  if(!capabilities("jpeg"))
  {
    return(false(no_capability_msg("jpeg")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_png_capability <- function()
{
  if(!capabilities("png"))
  {
    return(false(no_capability_msg("png")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_tiff_capability <- function()
{
  if(!capabilities("tiff"))
  {
    return(false(no_capability_msg("tiff")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_tcltk_capability <- function()
{
  if(!capabilities("tcltk"))
  {
    return(false(no_capability_msg("tcltk")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_x11_capability <- function()
{
  if(!capabilities("X11"))
  {
    return(false(no_capability_msg("X11")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_aqua_capability <- function()
{
  if(!capabilities("aqua"))
  {
    return(false(no_capability_msg("aqua")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_http_ftp_capability <- function()
{
  if(!capabilities("http/ftp"))
  {
    return(false(no_capability_msg("http/ftp")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_sockets_capability <- function()
{
  if(!capabilities("sockets"))
  {
    return(false(no_capability_msg("sockets")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_libxml_capability <- function()
{
  if(!capabilities("libxml"))
  {
    return(false(no_capability_msg("libxml")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_fifo_capability <- function()
{
  if(!capabilities("fifo"))
  {
    return(false(no_capability_msg("fifo")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_cledit_capability <- function()
{
  if(!capabilities("cledit"))
  {
    return(false(no_capability_msg("cledit")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_iconv_capability <- function()
{
  if(!capabilities("iconv"))
  {
    return(false(no_capability_msg("iconv")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_nls_capability <- function()
{
  if(!capabilities("NLS"))
  {
    return(false(no_capability_msg("NLS")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_profmem_capability <- function()
{
  if(!capabilities("profmem"))
  {
    return(false(no_capability_msg("profmem")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_cairo_capability <- function()
{
  if(!capabilities("cairo"))
  {
    return(false(no_capability_msg("cairo")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_icu_capability <- function()
{
  if(as.package_version(version) < "3.1.2")
  {
    return(
      false(not_declared_msg("ICU" , "3.1.2."))
    )
  }
  if(!capabilities("ICU"))
  {
    return(false(no_capability_msg("ICU")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_long_double_capability <- function()
{
  if(as.package_version(version) < "3.1.3")
  {
    return(
      false(not_declared_msg("long.double" , "3.1.3."))
    )
  }
  if(!capabilities("long.double"))
  {
    return(false(no_capability_msg("long.double")))
  }
  TRUE
}

#' @rdname r_has_jpeg_capability
#' @export
r_has_libcurl_capability <- function()
{
  if(as.package_version(version) < "3.2.0")
  {
    return(
      false(not_declared_msg("libcurl" , "3.2.0."))
    )
  }
  if(!capabilities("libcurl"))
  {
    return(false(no_capability_msg("libcurl")))
  }
  TRUE
}


no_capability_msg <- function(cap)
{
  gettextf("R does not have %s capability.", cap)
}

not_declared_msg <- function(cap, vrsn)
{
  gettextf("%s capability is not declared for versions of R before %s.", cap)
}