help.BRugs <- function(browser = getOption("browser"))
{
    ## stolen from help.start()
    #    a <- system.file("OpenBUGS", "Manuals", "WinBUGS Manual.html", package="BRugs")
    #    if (!file.exists(a))
    #        stop("I can't find the html help")
    #    a <- chartr("/", "\\", a)
    #    message("If nothing happens, you should open `", a, "' yourself")
    #    browseURL(a, browser = browser)
    #    invisible("")
    ## Andrew now omits the BRugs introduction, hence just pointing to help.WinBUGS these days:
    help.WinBUGS(browser = browser)
}

help.WinBUGS <- function(browser = getOption("browser"))
{
    # stolen from help.start()
    a <- file.path(options()$OpenBUGS, "Manuals", "Contents.html")
    if (!file.exists(a))
        stop("HTML help not found in file ", a)
    if (is.R())
      a <- chartr("/", "\\", a)
    else
      a <- gsub ("/", "\\\\", a)
    message("If nothing happens, you should open `", a, "' yourself")
    browseURL(a, browser = browser)
    invisible("")
}
