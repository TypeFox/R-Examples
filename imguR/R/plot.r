plot.imgur_image <-
plot.imgur_gallery_image <-
function(x, ...) {
    ext <- tolower(file_ext(x$link))
    if(ext == 'gif') {
        stop("plot method for GIF images not currently suported")
    } else {
        readfun <- switch(ext, 'jpg' = readJPEG, 
                               'jpeg' = readJPEG, 
                               'png' = readPNG,
                               stop('Unrecognized file extension.'))
        contents <- GET(x$link)
        out <- do.call(readfun, list(contents))
        plot(NULL, xlim=c(0,1), ylim=c(0,1), xaxt='n', yaxt='n', 
             xlab = "", ylab = "", main = "",
             xaxs='i', yaxs='i', mar=rep(0,4), mgp=rep(0,3), ...)
        rasterImage(out, 0,0,1,1)
    }
    invisible(x)
}
