tilePerim0 <- function (object,inclbdry=TRUE) {
    x <- object[["x"]]
    y <- object[["y"]]
    xx <- c(x,x[1])
    yy <- c(y,y[1])
    if(inclbdry) {
        ok <- rep(TRUE,length(x))
    } else {
        bp1 <- object[["bp"]]
        bp2 <- c(bp1,bp1[1])
        bpm <- cbind(bp1,bp2[-1])
        ok  <- !apply(bpm,1,all)
    }
    sum(sqrt(((xx[-1] - x)[ok])^2 + ((yy[-1] - y)[ok])^2))
}
