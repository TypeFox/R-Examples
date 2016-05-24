khr2estUDm <- function(x)
{
    if (!inherits(x, "khr"))
        stop("non convenient class")
    ii <- lapply(x, function(y) {
        ud <- asc2spixdf(y$UD)
        ud <- new("estUD", ud)
        if (y$hmeth == "href") {
            hli <- list(h = y$h, meth = "href")
        } else {
            if (y$hmeth == "LSCV") {
                hli <- y$h
                hli$meth <- "LSCV"
            } else {
                if (inherits(x, "khrud")|inherits(x, "khrvol")) {
                    hli <- list(h = y$h, meth = "specified")
                } else {
                    if (inherits(x, "kbbhrud")) {
                        hli <- list(values = y$h,
                                    meth = "BB-specified")
                    }
                }
            }
        }
        slot(ud, "h") <- hli
        if (inherits(x, "khrvol")) {
            slot(ud, "vol") <- TRUE
        } else {
            slot(ud, "vol") <- FALSE
        }
        names(slot(ud,"data"))[1] <- "ud"
        return(ud)
    })
    names(ii) <- names(x)
    class(ii) <- "estUDm"
    if (length(ii) < 2) {
        ii <- ii[[1]]
    }
    return(ii)
}


