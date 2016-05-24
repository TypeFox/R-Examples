plotnonp <- function(data, x=2, y=c(3,4,5,7,8), ylim=NULL, lty=c(1,2,3,2,3),
                     cols=rep(1,length(y)), month, year, step=12, xlab, ylab, ...)
{
    if(is.data.frame(data))
        block <- data
    else
        block <- read.table(data, header=TRUE)

    if(is.character(x))
        x <- match(x, dimnames(block)[[2]])
    if(is.character(y))
        y <- match(y, dimnames(block)[[2]])

    if(missing(xlab))
        xlab <- dimnames(block)[[2]][x]
    if(missing(ylab))
        ylab <- paste("f(", xlab, ")")

    if(is.null(ylim))
        ylim <- c(min(block[,y], na.rm=TRUE), max(block[,y], na.rm=TRUE))

    plot(block[,x], block[,y[1]], ylim=ylim, lty=lty[1], col=cols[1], type="l", 
         xlab=xlab, ylab=ylab, axes=FALSE, ...)
    box()
    axis(2)

    if(!missing(month) & !missing(year)){
        start <- block[1, x] - month + 1
        stop <- max(block[, x] + 1, na.rm=TRUE)
        pos <- seq(start, stop, step)
        label <- (pos - pos[1])/step + year
        if(nrow(block) <= 24) {
            if(step == 12) {
                label2 <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
            }
            else if(step == 4) {
                label2 <- c("Jan", "Apr", "Jul", "Oct")
            }
            else if(step == 2) {
                label2 <- c("Jan", "Jul")
            }
            else {
                label2 <- FALSE
            }
            label2 <- rep(label2, length.out = nrow(block) + month - 1)
            label2 <- label2[month:(nrow(block) + month - 1)]
            start2 <- block[1, x]
            stop2 <- max(block[, x], na.rm = TRUE)
            pos2 <- seq(start2, stop2, 1)
            axis(side = 1, at = pos2, labels = label2, mgp = c(3, 0.5, 0))
            axis(side = 1, at = pos, labels = label, mgp = c(3, 1.5, 0))
        }
        else axis(side = 1, at = pos, labels = label)
    }
    else axis(1)

    if(length(y) > 1){
        for(i in 2:length(y))
            lines(block[,x], block[,y[i]], lty=lty[i], col=cols[i])
    }

    return(invisible())   
}
