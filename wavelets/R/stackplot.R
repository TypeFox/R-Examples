stackplot <- function (x, y = NULL, type = "h", axes.labels = FALSE, xlab = "", ylab = "", y.llabs = NULL, y.rlabs = NULL, draw.divides = TRUE, xtick.at = NULL, ytick.at = NULL, col = "black", main = "")
{
    innerstackplot <- function(x, matrix, type, colorp, vertshift, plotarray)
    {
        plotmatrix <- matrix

        plot(x, plotmatrix[,1], type = type, xlim = xlimits, ylim = c(mindata[1], mindata[1] + totalrange), ylab = ylab,  xlab = xlab, axes = FALSE, frame.plot = TRUE, col = colorp)
    
        if(dim(plotmatrix)[2] > 1) {
            if(type == "h") {
                for (i in 2:dim(plotmatrix)[2])
                {
                    shifty1 <- vertshift[i]
                    shifty2 <- plotmatrix[,i] + vertshift[i]
                    segments(x, shifty1, x, shifty2, col = colorp)
                }
            }
            if(type == "l") {
                for (i in 2:dim(plotmatrix)[2])
                {
                    shifty2 <- plotmatrix[,i] + vertshift[i]
                    lines(x, shifty2, col = colorp)
                }
            }
        }
    }

    find.rightlabelpos <- function ()
    {
        right.usrplotrange <- par()$usr[2] - par()$usr[1]
        NDCplotrange <- par()$plt[2] - par()$plt[1] 
        marginpos <- (1-par()$plt[2])/2
        right.usrlabelpos <- ((marginpos*right.usrplotrange)/NDCplotrange) + par()$usr[2]
        right.usrlabelpos
    }

    draw.rightlabels <- function(y.rlabs, right.usrlabelpos, vertshift)
    {
        for (i in 1:length(y.rlabs)) 
        {
            text(right.usrlabelpos, vertshift[i], y.rlabs[[i]], xpd = TRUE)
        }
    }

    find.leftlabelpos <- function ()
    {
        left.usrplotrange <- par()$usr[1]
        NDCplotrange <- par()$plt[1]
        marginpos <- (par()$plt[1])/2
        left.usrlabelpos <- ((marginpos*left.usrplotrange)/NDCplotrange)
        left.usrlabelpos
    }

    draw.leftlabels <- function(y.llabs, left.usrlabelpos, vertshift)
    {
        for (i in 1:length(y.llabs)) 
        {
            text(left.usrlabelpos, vertshift[i], y.llabs[[i]], xpd = TRUE)
        }
    }

    coercetomatrix <- function(y) { 
        for(i in 1:length(y)) {
            if(class(y[[i]])[1] == "mts" || class(y[[i]])[1] == "ts") {

                x.start <- start(y[[i]])[1]
                x.end <- end(y[[i]])[1]
                x.range <- x.start:x.end
        
                if(class(y[[i]])[1] == "mts") {
                    y.matrix <- NULL
                    for(j in 1:length(y[[i]][1,])) {
                        y.matrix <- cbind(y.matrix, as.matrix(y[[i]][,j]))
                    }
                    y[[i]] <- y.matrix     
                }
                else {
                    y[[i]] <- as.matrix(y[[i]])
                }
                rownames(y[[i]]) <- x.range

            }
            if(class(y[[i]])[1] != "matrix") {
                y[[i]] <- as.matrix(y[[i]])
            }
        }

        y 
    }

    if(is.null(y) == TRUE) {
        y <- x
        if(class(y)[1] != "list") {
            y <- list(y)
        }

        y <- coercetomatrix(y)
           
        if(is.null(rownames(y[[1]])) == TRUE) {
             x <- list(1:length(y[[1]][,1]))
             if(length(y) != 1) {
                 for(k in 2:length(y)) {
                     x <- c(x, list(1:length(y[[k]][,1])))
                 }
             }
        }
        else { 
             x <- list(as.numeric(rownames(y[[1]])))
             if(length(y) != 1) {
                 for(k in 2:length(y))     {
                     if(is.null(rownames(y[[k]])) == TRUE) {
                         x <- c(x, list(1:length(y[[k]][,1])))
                     }
                     else {
                         x <- c(x, list(as.numeric(rownames(y[[k]]))))
                     } 
                 }
             }
        }
    }
    else {
        if(class(x) != "list") {
            x <- list(x)
        }
        if(class(y) != "list") {
            y <- list(y)
        }

        y <- coercetomatrix(y)

        if(length(x) != length(y)) {
            k <- 0
            add.type <- x[[1]]
            for(k in 1:(length(y) - length(x))) {
                x <- c(x, list(add.type))
            } 

        }

        for(k in 1:length(y)) {
            if(length(x[[k]]) != length(y[[k]][,1])) {
                stop("x and y lengths differ.")
            }
        }
    }

    if(length(x[[1]]) != length(y[[1]][,1])) {
        stop("x and y lengths differ.")
    }


    if(is.null(type[[1]] == TRUE)) {
        type[[1]] <- "h"
    }
    if(length(type) < length(y)) {
        for(i in 2:length(y)) { 
            type <- c(type, list(type[[i-1]]))
        }
    }

    if (length(col) < length(y)) {
        if(is.null(col[1] == TRUE)) {
            col[1] <- "black"
        }
        for(i in 2:length(y)) { 
            if(is.null(col[i]) == TRUE) {
                col <- c(col, col[i-1])
            }
        }
    }

    plotmatrix <- as.matrix(y[[1]])
    plotarray <- 1:dim(plotmatrix)[2]
  
    rangeperplot <- NULL 
    mindata <- NULL

    for(i in 1:length(plotmatrix[1,])) {
        if(min(plotmatrix[,i], na.rm = TRUE) > 0) { 
            range <- max(plotmatrix[,i], na.rm = TRUE)
            mindata <- c(mindata, 0)
            rangeperplot <- c(rangeperplot, range)     
        }
        if(min(plotmatrix[,i], na.rm = TRUE) < 0) {
            if(max(plotmatrix[,i], na.rm = TRUE) < 0) {
                range <- abs(min(plotmatrix[,i], na.rm = TRUE))
                mindata <- c(mindata, min(plotmatrix[,i], na.rm = TRUE))
                rangeperplot <- c(rangeperplot, range)
            }
            else {
                range <- max(plotmatrix[,i], na.rm = TRUE) - min(plotmatrix[,i], na.rm = TRUE)                
                mindata <- c(mindata, min(plotmatrix[,i], na.rm = TRUE))
                rangeperplot <- c(rangeperplot, range)
            }
        }
    }

    totalrange <- sum(rangeperplot)

    jstartline <- mindata[1]
    for (j in plotarray) {
        jstartline <- jstartline + rangeperplot[j]
        plotarray[j] <- jstartline
    }
   
    xlimits <- c(min(x[[1]]), max(x[[1]]))
    if(axes.labels == TRUE) {
        plot(0, 0, type = "n", xlim = xlimits, ylim = c(mindata[1], mindata[1] + totalrange),  ylab = ylab,  xlab = xlab, frame.plot = TRUE, main = main)     
    }
    else {
        plot(0, 0, type = "n", xlim = xlimits, ylim = c(mindata[1], mindata[1] + totalrange),  ylab = ylab,  xlab = xlab, axes = FALSE, frame.plot = TRUE, main = main)
    }

    if(is.null(xtick.at) == FALSE && is.null(xlab) == TRUE) {
        axis(1, at = axTicks(1, xtick.at), labels = FALSE)
    }
    if(is.null(xtick.at) == FALSE && is.null(xlab) == FALSE) {
        axis(1, at = axTicks(1, xtick.at), labels = FALSE)
    }
    if(is.null(xtick.at) == TRUE && is.null(xlab) == TRUE) {
        axis(1,  labels = FALSE) 
    }
    if(is.null(xtick.at) == TRUE && is.null(xlab) == FALSE) {
        axis(1, labels = FALSE)
    }

    if(is.null(ytick.at) == FALSE && is.null(ylab) == TRUE) {
        ytickrate <- (ytick.at[2] - ytick.at[1])/ytick.at[3]   
        interval <- (totalrange - mindata[1])/ytickrate
        axis(2, at = axTicks(2, c(mindata[1], totalrange, interval)), labels = FALSE)
    }
    if(is.null(ytick.at) == FALSE && is.null(ylab) == FALSE) {
        ytickrate <- (ytick.at[2] - ytick.at[1])/ytick.at[3]   
        interval <- (totalrange - mindata[1])/ytickrate
        axis(2, at = axTicks(2, c(mindata[1], totalrange, interval)), labels = FALSE)
    }
    if(is.null(ytick.at) == TRUE && is.null(ylab) == TRUE) {
        axis(2,  labels = FALSE) 
    }
    if(is.null(ytick.at) == TRUE && is.null(ylab) == FALSE) {
        axis(2, labels = FALSE)
    }

    vertshift = rep(0, dim(plotmatrix)[2])
    sum <- 0

    if(dim(plotmatrix)[2] > 1) {
        for (j in 2:dim(plotmatrix)[2])
        {
            if(max(plotmatrix[,(j-1)], na.rm = TRUE) < 0) { #so min is less than 0
                vertshift[j] <- sum + abs(min(plotmatrix[,j], na.rm=TRUE))                
                sum <- sum + abs(min(plotmatrix[,j], na.rm=TRUE))
            }
            if(min(plotmatrix[,(j)], na.rm = TRUE) > 0) { #so max is greater than 0
                vertshift[j] <- sum + max(plotmatrix[,(j-1)], na.rm=TRUE)                
                sum <- sum + max(plotmatrix[,(j-1)], na.rm=TRUE)
            }
            if(min(plotmatrix[,(j)], na.rm = TRUE) <= 0 && max(plotmatrix[,(j-1)], na.rm = TRUE) >= 0) {
                vertshift[j] <- sum + max(plotmatrix[,(j-1)], na.rm=TRUE) - min(plotmatrix[,j], na.rm=TRUE)
                sum <- sum + max(plotmatrix[,(j-1)], na.rm=TRUE) - min(plotmatrix[,j], na.rm=TRUE)
            } 
            if(draw.divides == TRUE) {
                abline(plotarray[j-1], 0)
            }  
        }
    }

    par(new = TRUE)
    for(i in 1:length(y)) {
        innerstackplot(x[[i]], y[[i]], type = type[[i]], colorp = col[i], vertshift, plotarray)
        par(new = TRUE)
    }    

    if(is.null(y.rlabs) == FALSE) {
        draw.rightlabels(y.rlabs,  find.rightlabelpos(), vertshift)
    }

    if(is.null(y.llabs) == FALSE) {
        draw.leftlabels(y.llabs, find.leftlabelpos(), vertshift)
    }

    par(new = FALSE)
}
