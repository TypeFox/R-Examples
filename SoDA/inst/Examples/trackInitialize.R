setMethod("initialize", "track",
    function(.Object, x, y, ...) {
        if(missing(y)) {
            y <- x; x <- seq(along=y)
        }
        callNextMethod(.Object, x = x, y = y, ...)
    })
