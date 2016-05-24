"image_plot_info" <-
function (...) 
{
    temp <- list(...)
    xlim <- NA
    ylim <- NA
    zlim <- NA
    if (is.list(temp[[1]])) {
        xlim <- range(temp[[1]]$x, na.rm = TRUE)
        ylim <- range(temp[[1]]$y, na.rm = TRUE)
        zlim <- range(temp[[1]]$z, na.rm = TRUE)
    }
    if (is.matrix(temp[[1]])) {
        xlim <- c(0, 1)
        ylim <- c(0, 1)
        zlim <- range(temp[[1]], na.rm = TRUE)
    }
    if (length(temp) >= 3) {
        if (is.matrix(temp[[3]])) {
            xlim <- range(temp[[1]], na.rm = TRUE)
            ylim <- range(temp[[2]], na.rm = TRUE)
            zlim <- range(temp[[3]], na.rm = TRUE)
        }
    }
    xthere <- match("x", names(temp))
    ythere <- match("y", names(temp))
    zthere <- match("z", names(temp))
    if (!is.na(zthere)) 
        zlim <- range(temp$z, na.rm = TRUE)
    if (!is.na(xthere)) 
        xlim <- range(temp$x, na.rm = TRUE)
    if (!is.na(ythere)) 
        ylim <- range(temp$y, na.rm = TRUE)
    if (!is.null(temp$zlim)) 
        zlim <- temp$zlim
    if (!is.null(temp$xlim)) 
        xlim <- temp$xlim
    if (!is.null(temp$ylim)) 
        ylim <- temp$ylim
    list(xlim = xlim, ylim = ylim, zlim = zlim)
}

