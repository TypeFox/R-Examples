filled.gap <-
function(dat, gcol='gray', edges = T, col = c('red','red'), 
                       lty = c(1,1), append = F, density = NULL, angle=45, ... ){
    dat <- as.data.frame(dat)
    # check the endings (no NAs)
    if (any(is.na(dat[1,]))) dat <- dat[-1,]
    if (any(is.na(dat[nrow(dat),]))) dat <- dat[-nrow(dat),]
    xx <- c(dat[[1]], rev(dat[[1]]))
    yy <- c(dat[[2]], rev(dat[[3]]))
    if (!append) plot(xx, yy, type='n', ...)
    polygon(xx, yy, border=NA, col=gcol, density=density, angle=angle)
    if (edges) matlines(dat[[1]], dat[2:3], col=col, lty=lty)
}
