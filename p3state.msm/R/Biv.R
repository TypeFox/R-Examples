`Biv` <-function (object, time1, time2) 
{
    if (missing(object)) 
        stop("Argument 'object' is missing with no default")
    if (inherits(object, "p3state")) 
        mydata <- object$datafr
    if (inherits(object, "data.frame")) 
        mydata <- object
    if (!inherits(object, "data.frame") & !inherits(object, "p3state")) 
        stop("'object' must be of class 'p3state'")
    if (missing(time1)) 
        stop("Argument 'time1' is missing with no default")
    if (missing(time2)) 
        stop("Argument 'time2' is missing with no default")
    if (time1 < 0 | time2 < 0) 
        stop("'time1' and 'time2' must be positive")
    p1 <- which(mydata[, 1] <= time1 & mydata[, 3] <= time2 & 
        mydata[, 2] * mydata[, 5] == 1)
    if (length(p1)==0) res<-0							#changed
    else res <- sum(mydata[p1, 6])						#changed
    return(res)
}