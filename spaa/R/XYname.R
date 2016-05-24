XYname <-
function(x,y){
    plotdata <- 1:(x*y)
    dim(plotdata) <- c(x,y)
    mat <- as.matrix(dist(plotdata))
    xxx <- expand.grid(1:x, 1:y)
    xlab <- paste("X",xxx[,1], sep = "")
    ylab <- paste("Y",xxx[,2], sep = "")
    plotlab <- paste(xlab, ylab, sep = "")
    return(plotlab)
}
