mypal <- function(set=TRUE,...) {
    oldpal <- palette()
    col <- c("black","darkblue","darkred","goldenrod","mediumpurple",
             "seagreen","aquamarine3","violetred1","salmon1",
             "lightgoldenrod1","darkorange2","firebrick1","violetred1", "gold")
    if (!set) return(col)
    palette(col)
    invisible(oldpal)
}


##' This function transforms a standard color (e.g. "red") into an
##' transparent RGB-color (i.e. alpha-blend<1).
##'
##' This only works for certain graphics devices (Cairo-X11 (x11 as of R>=2.7), quartz, pdf, ...).
##' @title Generate a transparent RGB color
##' @param col Color (numeric or character)
##' @param alpha Degree of transparency (0,1)
##' @param locate Choose colour (with mouse)
##' @return   A character vector with elements of 7 or 9 characters, '"\#"'
##'  followed by the red, blue, green and optionally alpha values in
##' hexadecimal (after rescaling to '0 ... 255').
##' @author Klaus K. Holst
##' @examples
##' plot(runif(1000),cex=runif(1000,0,4),col=Col(c("darkblue","orange"),0.5),pch=16)
##' @keywords color
##' @export
Col <- function(col,alpha=0.2,locate=0) {
    if (locate>0) colsel(locate)
    
    mapply(function(x,alpha)
        do.call(rgb,as.list(c(col2rgb(x)/255,alpha))),
        col,alpha)
}

