.onLoad <- function (lib, pkg)
{
    ## Starting the DDE server automatically if under Windows
    ## and option use.DDE == TRUE
    use.DDE <- getOption("use.DDE")
    if (.Platform$OS.type == "windows" && !is.null(use.DDE) && use.DDE)
        guiDDEInstall()

    ## If an IDE is defined, start it now
    IDE <- getOption("IDE")
    if (!is.null(IDE) && file.exists(IDE))
        system(paste("\"", IDE, "\"", sep = ""), wait = FALSE)
}

### This is to define writeClipboard() on other platforms than windows (but does
### nothing currently!)
if (.Platform$OS.type != "windows")
    "writeClipboard" <- function (str, format = 1)
        stop("Not implemented yet on other platforms than Windows")
