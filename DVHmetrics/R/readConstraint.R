readConstraint <-
function(x, ...) {
    fn <- if(!missing(x)) {
        Sys.glob(x)
    } else if(interactive()) {
        file.choose()
    } else {
        NULL
    }

    if(!is.null(fn)) {
        x <- read.table(file=fn, header=TRUE, stringsAsFactors=FALSE, ...)
        cNames <- tolower(names(x))
        cNames[cNames == "patid"] <- "patID"
        setNames(x, cNames)
    } else {
        warning("No single file was selected")
        NULL
    }
}
