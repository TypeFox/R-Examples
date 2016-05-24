 image_plot_plt <-
function (...,add = FALSE, legend.shrink = 0.9, legend.width = 0.04, 
    horizontal = FALSE, offset = 2 * legend.width, bigplot = NULL, 
    smallplot = NULL) 
{
    old.par <- par(no.readonly = TRUE)
    if (is.null(smallplot)) 
        stick <- TRUE
    else stick <- FALSE
    if (is.null(smallplot)) {
        smallplot <- old.par$plt
        if (horizontal) {
            smallplot[3] <- offset
            smallplot[4] <- smallplot[3] + legend.width
            pr <- (smallplot[2] - smallplot[1]) * ((1 - legend.shrink)/2)
            smallplot[1] <- smallplot[1] + pr
            smallplot[2] <- smallplot[2] - pr
        }
        else {
            smallplot[2] <- 1 - offset
            smallplot[1] <- smallplot[2] - legend.width
            pr <- (smallplot[4] - smallplot[3]) * ((1 - legend.shrink)/2)
            smallplot[4] <- smallplot[4] - pr
            smallplot[3] <- smallplot[3] + pr
        }
    }
    if (is.null(bigplot)) {
        bigplot <- old.par$plt
        if (!horizontal) {
            bigplot[2] <- min(bigplot[2], smallplot[1] - offset)
        }
        else {
            bottom.lines <- old.par$mar[1]
            if (bottom.lines > 0) {
                off2 <- (3 * bigplot[1])/bottom.lines
            }
            else {
                off2 <- 0
            }
            bigplot[3] <- max(bigplot[3] - off2, smallplot[4] + 
                offset)
        }
    }
    if (stick) {
        if (!horizontal) {
            dp <- smallplot[2] - smallplot[1]
            smallplot[1] <- min(bigplot[2] + offset, smallplot[1])
            smallplot[2] <- smallplot[1] + dp
        }
    }
    return(list(smallplot = smallplot, bigplot = bigplot))
}
