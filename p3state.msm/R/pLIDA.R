pLIDA <- function (object, time1, time2, tp = NULL) {
    if (missing(object))
        stop("Argument 'object' is missing with no default")
    if (inherits(object, "p3state"))
        mydata <- object$datafr
    if (inherits(object, "data.frame"))
        mydata <- object
    if (!inherits(object, "data.frame") & !inherits(object, "p3state"))
        stop("'object' must be of class 'p3state'")
    if (missing(time1))
        time1 <- 0
    if (missing(tp))
        tp <- "all"
    if (missing(time2))
        stop("Argument 'time2' is missing with no default")
    if (time1 < 0 | time2 < 0 | time1 > time2)
        stop("'time1' and 'time2' must be positive, and time1 < time2")
    if (tp == "all" | tp == "p11") {
        p1 <- max(which(mydata[, 1] <= time2))
        p2 <- max(which(mydata[, 1] <= time1))
        aux1 <- mydata[p1, 7]
        aux2 <- mydata[p2, 7]
        if (aux2 == 0)
            restp11 <- 0
        else restp11 <- aux1/aux2
    }
    if (tp == "all" | tp == "p12") {
	p01 <- max(which(mydata[, 1] <= time2))
	p1 <- which(mydata[, 1] <= time2 & mydata[, 1] > time1 &
            mydata[, 4] <= time2 & mydata[, 5] == 1)
        p2 <- max(which(mydata[, 1] <= time1))
	aux01 <- mydata[p01, 7]
	aux2 <- mydata[p2, 7]
        if (aux2 == 0)
            restp12 <- 0
        else restp12 <- (aux2-aux01-sum(mydata[p1, 6]))/(aux2)
        if (restp12 < 0) restp12 <- 0
        if (restp12 > 1) restp12 <- 1
    }
    if (tp == "all" | tp == "p22") {
	p01 <- max(which(mydata[, 1] <= time1))
	p1 <- which(mydata[, 1] <= time1 & mydata[, 4] <= time2 &
            mydata[, 5] == 1)
        p2 <- which(mydata[, 1] <= time1 & mydata[, 4] <= time1 &
            mydata[, 5] == 1)
	aux1 <- mydata[p01, 7]
        restp22 <- (1-aux1-sum(mydata[p1, 6]))/(1-aux1-sum(mydata[p2, 6]))
        if (restp22 < 0) restp22 <- 0
        if (restp22 > 1) restp22 <- 1
    }
    if (tp == "p11")
        res <- restp11
    if (tp == "p12")
        res <- restp12
    if (tp == "p22")
        res <- restp22
    if (tp == "all")
        res <- list(restp11, restp12, restp22)
    return(res)
}
