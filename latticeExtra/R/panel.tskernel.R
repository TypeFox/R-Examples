
panel.tskernel <-
    function(x, y, ..., 
             width = NROW(x) %/% 10 + 1, n = 300,
             c = 1, sides = 2, circular = FALSE,
             kern = kernel("daniell", rep(floor((width/sides)/sqrt(c)), c)))
{
    if (!missing(kern))
        .Deprecated("The 'kern' argument to panel.tskernel is Deprecated. Use simpleSmoothTs directly.")
    if (!missing(y)) {
        x <- as.numeric(x)
        y <- as.numeric(y)
        stopifnot(NCOL(x) == 1)
        if (diff(range(diff(x))) > getOption("ts.eps"))
            stop("'x' should be a regular series")
        x <- ts(y, start = x[1], end = tail(x,1), deltat = diff(x[1:2]))
    }
    x <- as.ts(x)
    s <- simpleSmoothTs(x, width = width, c = c, sides = sides,
                        circular = circular, kern = kern, n = n)
    panel.lines(s, ...)
}
