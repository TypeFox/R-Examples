run <- function(n = 100, ...) {
    # direct HTML output
    out("<H2>Some HTML</h2>")
    # all arguments are passed as strings from the URL, so convert to numeric as needed
    n <- as.integer(n)
    # create a WebPlot
    p <- WebPlot(800, 600)
    x <- runif(n)
    plot(x, rnorm(n), pch=19, col=2)
    # insert the plot in the page
    out(p)
    # verbatim print
    oprint(n)
    oprint(summary(x))
    # HTML table
    otable(data.frame(a=1:5, b=c("a","b","c","d","e")))
    # return the whole page
    done()
}
