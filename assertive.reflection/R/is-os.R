#' @rdname is_windows
#' @export
is_bsd <- function()
{
  if(!grepl("BSD", Sys.info()[["sysname"]]))
  {
    return(not_this_os("BSD-based"))
  }
  TRUE
}

#' @rdname is_windows
#' @export
is_linux <- function()
{
  if(Sys.info()["sysname"] != "Linux")
  {
    return(not_this_os("Linux"))
  }
  TRUE
}

#' @rdname is_windows
#' @export
is_mac <- function()
{
  if(Sys.info()["sysname"] != "Darwin")
  {
    return(not_this_os("OS X"))
  }
  TRUE
}

#' @rdname is_windows
#' @export
is_osx <- is_mac
  
#' @rdname is_windows
#' @export
is_osx_cheetah <- function()
{
  is_osx_version("Cheetah")
}

#' @rdname is_windows
#' @export
is_osx_puma <- function()
{
  is_osx_version("Puma")
}

#' @rdname is_windows
#' @export
is_osx_jaguar <- function()
{
  is_osx_version("Jaguar")
}

#' @rdname is_windows
#' @export
is_osx_panther <- function()
{
  is_osx_version("Panther")
}

#' @rdname is_windows
#' @export
is_osx_tiger <- function()
{
  is_osx_version("Tiger")
}

#' @rdname is_windows
#' @export
is_osx_leopard <- function()
{
  is_osx_version("Leopard")
}

#' @rdname is_windows
#' @export
is_osx_snow_leopard <- function()
{
  is_osx_version("Snow Leopard")
}

#' @rdname is_windows
#' @export
is_osx_lion <- function()
{
  is_osx_version("Lion")
}

#' @rdname is_windows
#' @export
is_osx_mountain_lion <- function()
{
  is_osx_version("Mountain Lion")
}

#' @rdname is_windows
#' @export
is_osx_mavericks <- function()
{
  is_osx_version("Mavericks")
}

#' @rdname is_windows
#' @export
is_osx_yosemite <- function()
{
  is_osx_version("Yosemite")
}

#' @rdname is_windows
#' @export
is_osx_el_capitan <- function()
{
  is_osx_version("El Capitan")
}

#' @rdname is_windows
#' @export
is_solaris <- function()
{
  if(Sys.info()["sysname"] != "SunOS")
  {
    return(not_this_os("Solaris"))
  }
  TRUE
}

#' @rdname is_windows
#' @export
is_unix <- function()
{
  if(.Platform$OS.type != "unix")
  {
    return(not_this_os("Unix-based"))
  }
  TRUE
}

#' What OS is running?
#' 
#' Is the operating system in this machine Windows/Unix/Mac based.
#' 
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_windows} returns \code{TRUE} if the OS on the current 
#' platform is Microsoft windows-based.  \code{is_unix} returns \code{TRUE} if 
#' the OS is Unix based (pretty much anything that isn't Windows, including OS 
#' X). 
#' \code{is_mac}, \code{is_linux}, \code{is_bsd}, \code{is_solaris} return 
#' \code{TRUE} if the OS is Apple OS X, Linux, FreeBSD/NetBSD, or Solaris 
#' respectively.
#' \code{is_64_bit_os} returns \code{TRUE} when the operating system is 64-bit.
#' The \code{assert_*} functions return nothing but throw an error if the 
#' corresponding \code{is_*} functions return \code{FALSE}.
#' @references With the exception of \code{is_windows} and \code{is_unix} that 
#' use \code{.Platform$OS.type}, the OS is determined from 
#' \code{Sys.info()[["sysname"]]}, which (not on Windows) is calculated via the 
#' OS \code{uname} program.  GNU has more information on the return value: 
#' \url{https://www.gnu.org/software/libc/manual/html_node/Platform-Type.html}
#' and Wikipedia has a nice list of possible values: 
#' \url{https://en.wikipedia.org/wiki/Uname#Examples}
#' The names for different versions of Windows are decribed in:
#' \url{http://svn.r-project.org/R/trunk/src/library/utils/src/windows/util.c}
#' @seealso \code{\link[base]{.Platform}}, \code{\link[base]{Sys.info}}, 
#' \code{\link[base]{version}}, and \code{win.version}.
#' @examples
#' is_unix()
#' is_linux()
#' is_bsd()
#' is_solaris()
#' if(is_windows())
#' {
#'   is_windows_vista()
#'   is_windows_7()
#'   is_windows_8()
#'   is_windows_8.1()
#'   is_windows_10()
#'   is_windows_server_2008()
#'   is_windows_server_2008_r2()
#'   is_windows_server_2012()
#'   is_windows_server_2012_r2()
#' }
#' if(is_osx()) # is_mac is a synonym
#' {
#'   is_osx_cheetah()
#'   is_osx_puma()
#'   is_osx_jaguar()
#'   is_osx_panther()
#'   is_osx_tiger()
#'   is_osx_leopard()
#'   is_osx_snow_leopard()
#'   is_osx_lion()
#'   is_osx_mountain_lion()
#'   is_osx_mavericks()
#'   is_osx_yosemite()
#'   is_osx_el_capitan()
#' }
#' is_32_bit()
#' is_64_bit()
#' assertive.base::dont_stop(assert_is_windows())
#' assertive.base::dont_stop(assert_is_unix())
#' @export
is_windows <- function()
{
  if(.Platform$OS.type != "windows")
  {
    return(not_this_os("Windows"))
  }
  TRUE
}

#' @rdname is_windows
#' @export
is_windows_vista <- function()
{
  is_windows_version("Vista")
}

#' @rdname is_windows
#' @export
is_windows_7 <- function()
{
  is_windows_version("7")
}

#' @rdname is_windows
#' @export
is_windows_8 <- function()
{
  is_windows_version(">= 8")
}

#' @rdname is_windows
#' @export
is_windows_8.1 <- function()
{
  is_windows_version("8.1")
}

#' @rdname is_windows
#' @export
is_windows_10 <- function()
{
  is_windows_version("10")
}

#' @rdname is_windows
#' @export
is_windows_server_2008 <- function()
{
  is_windows_version("Server 2008")
}

#' @rdname is_windows
#' @export
is_windows_server_2008_r2 <- function()
{
  is_windows_version("Server 2008 R2")
}

#' @rdname is_windows
#' @export
is_windows_server_2012 <- function()
{
  is_windows_version("Server >= 2012")
}

#' @rdname is_windows
#' @export
is_windows_server_2012_r2 <- function()
{
  is_windows_version("Server 2012 R2")
}

#' Failure for bad OS
#' 
#' Wrapper to \code{false} for failure messages when the OS is not as 
#' expected.
#' @param os A string giving the name of the OS that was desired.
#' @return A string showing the results of \code{.Platform$OS} and 
#' \code{Sys.info()['sysname']}.
#' @seealso \code{\link[base]{.Platform}} and \code{\link[base]{Sys.info}}
#' @examples
#' \donttest{
#' assertive.reflection:::not_this_os("Windows")
#' assertive.reflection:::not_this_os("BSD-based")
#' }
#' @noRd
not_this_os <- function(os)
{
  false(
    gettext(
      "The operating system is not %s. R reports it as: Sys.info()['sysname'] = %s, .Platform$OS = %s."
    ), 
    os, 
    Sys.info()["sysname"],
    .Platform$OS
  )
}

is_windows_version <- function(version)
{
  if(!(ok <- is_windows()))
  {
    return(ok)
  }
  windows_name_text <- utils::win.version()
  windows_version <- sub(
    "Windows (10|Vista|7|>= 8|8.1|Server 2008|Server 2008 R2|Server >= 2012|Server >= 2012 R2).*", 
    "\\1",
    windows_name_text
  )
  if(windows_version != version)
  {
    return(
      false(
        gettext("The operating system is not Windows %s. R reports it as: Windows %s."), 
        version,
        windows_version
      )
    )
  }
  TRUE
}

is_osx_version <- function(version)
{
  if(!(ok <- is_mac()))
  {
    return(ok)
  }
  mac_version_text <- system("sw_vers -productVersion", intern = TRUE)
  mac_version <- as.numeric_version(mac_version_text)
  minor_version <- unlist(mac_version)[2]
  os_name <- switch(
    as.character(minor_version),
    "0" = "Cheetah",
    "1" = "Puma",
    "2" = "Jaguar",
    "3" = "Panther",
    "4" = "Tiger",
    "5" = "Leopard",
    "6" = "Snow Leopard",
    "7" = "Lion",
    "8" = "Mountain Lion",
    "9" = "Mavericks",
    "10" = "Yosemite",
    "11" = "El Capitan"
  )
  if(version != os_name)
  {
    return(
      false(
        gettext("The operating system is not OS X %s. R reports it as: OS X %s."), 
        version,
        os_name
      )
    )
  }
  TRUE
}
