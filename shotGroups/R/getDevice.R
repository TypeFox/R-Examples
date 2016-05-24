## open drawing device depending on context
getDevice <-
function() {
    ## if under RStudio -> don't open new device
    dev <- getOption("device")
    if(class(dev) == "character") {
        if(dev == "RStudioGD") {
            dev <- function(...) {}
            ## alternative: open new device explicitly depending on platform
            ## find out which function to use for opening a new plot window
            ## find out operating system
#             osType <- .Platform$OS.type
#             if(osType == "windows") {
#                 dev <- windows
#             } else if(osType == "unix") {
#                 sysName <- Sys.info()[["sysname"]]
#                 if(sysName == "Linux") {
#                     dev <- x11
#                 } else if(sysName == "Darwin") {
#                     dev <- quartz
#                 }
#             }
        }
    }                                    # if(class(dev) == "character")

    ## if not in interactive mode -> don't open new device
    if(!interactive()) { dev <- function(...) {} }
    return(dev)
}
