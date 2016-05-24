user <- function() as.vector(Sys.info()["user"])

runningWindows <- function() .Platform$OS.type == "windows"
runningLinux   <- function() Sys.info()["sysname"] == "Linux"

