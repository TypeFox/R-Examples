`rgl.isomap` <-
    function(x, web = "white", ...)
{
    if (!inherits(x, "isomap"))
        stop("'x' must be an 'isomap' result object")
    ordirgl(x, ...)
    z <- scores(x, ...)
    net <- x$net
    ## skip if web = NA
    if (any(!is.na(web))) {
        ## web can be a vector for points (or not): recycle
        web <- rep(web, length = nrow(z))
        if (is.factor(web))
            web <- as.numeric(web)
        web <- col2rgb(web)/255
        for (i in 1:nrow(net)) {
            lcol <- rgb(t(rowMeans(web[, net[i,]])))
            rgl.lines(z[net[i,],1], z[net[i,],2], z[net[i,],3], color=lcol, ...)
        }
    }
}
