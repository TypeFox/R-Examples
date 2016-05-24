lpr <- function (object, file = "Rplotlpr.ps", ...) 
{
    if (missing(object)) {
        current.device <- dev.cur()
        dev.off(dev.copy(device = postscript, file = file, ...))
        dev.set(current.device)
        system(paste("lpr", file))
        print(paste(file, "printed."))
    }
    else {
        if (missing(file)) 
            file <- "Robjlpr.txt"
        sink(file)
        object <- as.character(substitute(object))
        print(get(object))
        sink()
        system(paste("lpr", file))
        print(paste(object, "printed."))
    }
}
