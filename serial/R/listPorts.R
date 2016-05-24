#' Lists the serial interfaces.
#' 
#' This function lists all installed serial interfaces in a computer.
#' Thereby Windows, Linux and MacOs behave a little bit different. Please ensure
#' that you have the appropriate permissions to do a search in the registry or
#' in the corresponding linux folders.
#' 
#' @usage listPorts()
#' 
#' @section Windows:
#' 
#' In a Windows environment, this function tries to read out the registry keys
#' located in:
#' 
#' \code{"HKEY_LOCAL_MACHINE\\HARDWARE\\DEVICEMAP\\SERIALCOMM"}
#' 
#' This should be consistent with all installed hardware ports plus all virtual
#' ports.
#' 
#' @section Linux and MacOS:
#' 
#' Here the situation is a bit different, compared to Windows. All possible serial 
#' devices are located in \code{"/dev/tty[...]"} as a file connection. Still, all
#' virtual and closed dev's can be found here. This is confusing, because one will 
#' find more devices in this folder than physically (virtual) present.
#' In addition to that, on Ubuntu linux systems in "\code{/sys/devices/pnp0/...}"
#' only the plug and play devices of interest are listed again. That is the reason why, 
#' the function returns a subset of \code{"/dev/tty[...]"}, which is also present in the 
#' \code{"../pnp0/.."} folder.
#' 
#' On MacOs the installed interfaces are marked by "\code{tty.<name>}", with a unique 
#' name after the dot, which makes it easier to search for installed devices.
#' 
#' Subsequently, the user must know which interface is present and which isn't. AND the user
#' must have at least reading permissions in the corresponding folders. So in the end,
#' this function is a best guess of what is installed.
#' 
#' @return A character vector with the list of comports is returned.
#' @importFrom utils tail 
#' @export
listPorts <- function()
{
  sList <- ""
  .Tcl("package require platform")
  # get Platform, especially to ask for macOS
  tk_platform <- tclvalue(.Tcl("lindex [split [platform::generic] -] 0"))
  
  if(.Platform$OS.type == "windows")
  {
    .Tcl("package require registry")
    regpath <- paste("HKEY_LOCAL_MACHINE","HARDWARE","DEVICEMAP","SERIALCOMM",sep="\\\\")
    ser_devs <- tclvalue( .Tcl( paste("registry values",regpath) ) )
    ser_devs <- strsplit( ser_devs, " ")[[1]]
    sList <- sapply(ser_devs,function(val) tclvalue(.Tcl(paste("registry get",regpath,val))))
    attr(sList,"names")<-NULL
  }
  
  if(.Platform$OS.type == "unix") # the R's platform is more general
  {
    # get all possible tty's
    sList <- dir("/dev/",pattern = "tty[0SU.'ACM''USB']")
    
    if(tk_platform != "macosx") # any other unix
    {
      ## try to figure out the available ports
      w <- options(warn=-1)$warn # disable warnings
      p <- system("find /sys/ -iname tty\\* 2>/dev/null",intern=T)
      options(warn = w) # enable old sate
      
      # only Plug&Play dev's
      p <- p[grep("pnp",p)]
      # get the last entry of the string split
      p <- lapply(strsplit(p,"/"),tail,1)
      if(length(p[[1]]) == 0)
        p <- 0
      else
        p <- unlist(p)
      p <- p[!duplicated(p)]
      
      # only, if there is a match, in sList and 
      if(sum(p %in% sList) > 0)
        sList <- sList[sList %in% p]
    }
    message("Hint: On unix-systems the port list might be incomplete!")
  }
  
  if( length(sList) == 0 )
    sList <- ""
  return(sList[order(sList)])
}
