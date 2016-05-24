plot.p3state<-function (x, plot.trans = NULL, plot.marginal = NULL, plot.bivariate = NULL,
 time1, time2, xlab, ylab, zlab, col, col.biv = NULL,...) 
{
    if (missing(x)) 
        stop("Argument 'x' is missing with no default")
    if (!inherits(x, "p3state")) 
        stop("'x' must be of class 'p3state'")
     
    mydata <- x$datafr
    xlab.aux <- FALSE
    ylab.aux <- FALSE
    pmar <- FALSE
    pbiv <- FALSE
    if (missing(xlab)) {
        xlab <- "Time"
        xlab.aux <- TRUE
    }
    if (missing(ylab)) {
        ylab.aux <- TRUE
    }
    if (missing(plot.trans)) 
        plot.trans <- FALSE
    if (missing(plot.marginal)) 
        plot.marginal <- FALSE
    if (missing(plot.bivariate)) 
        plot.bivariate <- FALSE
    if (missing(col.biv)) 
        col.biv <- FALSE
    if (missing(time1)) 
        time1 <- 0
    if (missing(time2)) 
        time2 <- max(mydata[, 1])
    if (time1 < 0 | time2 < 0) 
        stop("'time1' and 'time2' must be positive")
    if (time1 > time2) 
        stop("Argument 'time1' cannot be greater then 'time2'")
    if (sum(mydata[, 2]) == 0) 
        stop("This is not a multi-state model")
    if (sum(mydata[, 2] == 0 & mydata[, 5] == 1) > 0) 
        ntrans <- 3
    else ntrans <- 2
    if (plot.bivariate == TRUE) 
        pbiv <- TRUE
    if (plot.marginal == TRUE) 
        pmar <- TRUE
    if (ntrans == 3) {
        plot.marginal <- FALSE
        plot.bivariate <- FALSE
    }
    if (plot.trans == "P11" | plot.trans == "all") {
        x <- which((mydata[, 2] == 1 | mydata[, 5] == 1) & mydata[, 
            1] > time1)
        y1 <- sort(mydata[x, 1])
        y2 <- unique(y1)
        y <- c(time1, y2)
        mydata2 <- matrix(data = 0, ncol = 2, nrow = length(y))
        for (k in 1:length(y)) {
            mydata2[k, 1] <- y[k]
        }
        for (k in 1:length(y)) {
            mydata2[k, 2] <- pLIDA(mydata, time1, mydata2[k, 
                1], tp = "p11")
        }
        if (ylab.aux == TRUE) 
            ylab <- paste("Estimated prob. p11(", time1, ",t)")
        matplot(mydata2[, 1], mydata2[, 2], xlab = xlab, ylab = ylab, 
            type = "s", ...)
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    if (plot.trans == "P12" | plot.trans == "all") {
        x <- which(mydata[, 1] > time1 & mydata[, 5] == 1)
        y1 <- c(mydata[x, 1] + 1e-09, mydata[x, 4] + 1e-09)
        y2 <- unique(y1)
        y <- sort(y2)
        y <- c(time1, y)
        mydata2 <- matrix(data = 0, ncol = 2, nrow = length(y))
        for (k in 1:length(y)) {
            mydata2[k, 1] <- y[k]
        }
        for (k in 1:length(y)) {
            mydata2[k, 2] <- pLIDA(mydata, time1, mydata2[k, 
                1], tp = "p12")
        }
        if (ylab.aux == TRUE) 
            ylab <- paste("Estimated prob. p12(", time1, ",t)")
        matplot(mydata2[, 1], mydata2[, 2], xlab = xlab, ylab = ylab, 
            type = "s", ...)
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    if (plot.trans == "P22" | plot.trans == "P23" | plot.trans == 
        "all") {
        x <- which(mydata[, 5] == 1 & mydata[, 1] <= time1)
        q2 <- c(mydata[x, 1], mydata[x, 4])
        q3 <- q2[q2 >= time1]
        y1 <- sort(q3)
        y2 <- unique(y1)
        y <- c(time1, y2)
        mydata2 <- matrix(data = 0, ncol = 2, nrow = length(y))
        for (k in 1:length(y)) {
            mydata2[k, 1] <- y[k]
        }
        for (k in 1:length(y)) {
            mydata2[k, 2] <- pLIDA(mydata, time1, mydata2[k, 
                1], tp = "p22")
        }
        if (ylab.aux == TRUE) 
            ylab <- paste("Estimated prob. p22(", time1, ",t)")
        if (plot.trans == "P22" | plot.trans == "all") 
            matplot(mydata2[, 1], mydata2[, 2], xlab = xlab, 
                ylab = ylab, type = "s", ...)
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
        if (ylab.aux == TRUE) 
            ylab <- paste("Estimated prob. p23(", time1, ",t)")
        if (plot.trans == "P23" | plot.trans == "all") 
            matplot(mydata2[, 1], 1 - mydata2[, 2], xlab = xlab, 
                ylab = ylab, type = "s", ...)
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    if (plot.marginal == "TRUE") {
        x <- which(mydata[, 5] == 1)
        y1 <- mydata[x, 4]
        y2 <- unique(y1)
        y3 <- sort(y2)
        y <- c(0, y3)
        mydata2 <- matrix(data = 0, ncol = 2, nrow = length(y))
        for (k in 1:length(y)) {
            mydata2[k, 1] <- y[k]
        }
        for (k in 1:length(y)) {
            mydata2[k, 2] <- Biv(mydata, max(mydata[, 1]), mydata2[k, 
                1])
        }
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
        if (ylab.aux == TRUE) 
            ylab <- "Estimated marginal dist. F2(t)"
        matplot(mydata2[, 1], mydata2[, 2], xlab = xlab, ylab = ylab, 
            type = "s", ...)
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    if (plot.bivariate == "TRUE") {
        x <- which(mydata[, 5] == 1)
        x1 <- unique(sort(mydata[x, 1]))
        y1 <- unique(sort(mydata[x, 3]))
        mydata2 <- matrix(data = 0, nrow = length(x1), ncol = length(y1))
        for (k in 1:length(x1)) {
            for (j in 1:length(y1)) mydata2[k, j] <- Biv(mydata, 
                x1[k], y1[j])
        }
        z <- mydata2
        if (xlab.aux == TRUE) 
            xlab <- "Time in state 1"
        if (ylab.aux == TRUE) 
            ylab <- "Time in state 2"
        if (missing(zlab)) 
            zlab <- "F(t1,t2)"
        if (missing(col)) 
            col <- "lightblue"
        persp(x1, y1, z, theta = 30, phi = 30, ticktype = "detailed", 
            expand = 0.5, xlab = xlab, ylab = ylab, zlab = zlab, 
            col = col, shade = 0.2, ...)
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
        x1 <- seq(0, max(mydata[mydata[, 5] == 1, 1]) * 1.25, 
            length.out = 30)
        y1 <- seq(0, max(mydata[mydata[, 5] == 1, 3]) * 1.25, 
            length.out = 30)
        data <- expand.grid(x1, y1)
        z <- seq(0, 1, length.out = 900)
        for (k in 1:900) z[k] <- Biv(mydata, data[k, 1], data[k, 
            2])
        z1 <- matrix(z, 30)
        
        if (col.biv == FALSE) {
        bw <- colours()[350-3*0:19]
        filled.contour(x1, y1, z1, xlab = "Time in state 1", 
            ylab = "Time in state 2",col = bw, ...)   }
              
        else {
        filled.contour(x1, y1, z1, xlab = "Time in state 1", 
            ylab = "Time in state 2", ...) }   
            
    }
    if (ntrans == 3 & pmar == TRUE) 
        cat("The plot for the marginal distribution of the second time cannot be given for the illness-death model", 
            "\n")
    if (ntrans == 3 & pbiv == TRUE) 
        cat("The plot for the bivariate distribution function cannot be given for the illness-death model", 
            "\n")
}



