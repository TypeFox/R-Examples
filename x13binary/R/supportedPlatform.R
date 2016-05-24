#' Test Platform Support
#' 
#' Returns \code{TRUE} is the platform is supported, \code{FALSE} otherwise.
#'  
#' @examples
#' supportedPlatform()
#' 
#' @export
#' @import utils
supportedPlatform <- function(){

  z <- FALSE  # in case we missed something

  # Windows
  if (.Platform$OS.type == "windows"){
    z <- TRUE
  }

  # Linux
  if (Sys.info()["sysname"] %in% c("Linux")){
    z <- TRUE
  }

  # OS-X
  if (Sys.info()["sysname"] %in% c("Darwin")){
    # Darwin Version numbers for OS-X Versions 
    # https://en.wikipedia.org/wiki/Darwin_(operating_system)
    #     10  OS X Snow Leopard
    # 11.0.0  OS X Lion
    # 12.0.0  OS X Mountain Lion
    # 12.6.0  
    # 13.0.0  OS X Mavericks
    # 13.4.0  
    # 14.0.0  OS X Yosemite
    # 14.5.0  
    # 15.0.0  OS X El Capitan
    # 15.2.0  

    # OS-X needs to be at least Lion (TODO try out)
    z <- compareVersion(Sys.info()["release"], "11.0.0") >= 0
  }

  # Other Unix (eg Solaris)
  if ((.Platform$OS.type == "unix") && 
      !(Sys.info()["sysname"] %in% c("Darwin", "Linux"))){
    z <- FALSE
  }

  z
}
