partana <- function (c, dist) 
{
   UseMethod("partana")
}

partana.default <- function (c, dist) 
{
    c <- as.integer(clustify(c))
    numclu <- max(c)
    call <- match.call()
    if (class(dist) != 'dist') {
        stop('you must supply an object of class dist')
    }
    x <- (max(1,max(dist)) - as.matrix(dist)) / max(1,max(dist))
    num <- 0
    sumnum <- 0
    den <- 0
    sumden <- 0
    card <- rep(0,numclu)
    for (i in 1:numclu) {
        card[i] <- sum(c==i)
    }
    numplt <- nrow(x)
    simptc <- matrix(0, nrow = numplt, ncol = numclu)
    simctc <- matrix(0, nrow = numclu, ncol = numclu)
    for (i in 1:numplt) {
        for (j in 1:numclu) {
            if (c[i] == j) {
                if (card[j] > 1) {
                  simptc[i, j] <- (sum(x[i, c == j]) - 1)/(card[j] - 
                    1)
                }
                else {
                  simptc[i, j] <- 1
                }
            }
            else {
                if (card[j] > 0) {
                  simptc[i, j] <- sum(x[i, c == j])/card[j]
                }
                else {
                  simptc[i, j] <- 0
                }
            }
        }
    }
    for (i in 1:numclu) {
        for (j in 1:numclu) {
            if (i == j) {
                if (card[i] > 1) {
                  simctc[i, j] <- (sum(x[c == i, c == i]) - card[i])/(card[i]^2 - 
                    card[i])
                  sumnum <- sumnum + simctc[i, i] * ((card[i]^2 - 
                    card[i])/2)
                  num <- num + ((card[i]^2 - card[i])/2)
                }
            }
            else {
                if (card[i] != 0 & card[j] != 0) {
                  simctc[i, j] <- sum(x[c == i, c == j])/(card[i] * 
                    card[j])
                  sumden <- sumden + sum(x[c == i, c == j])
                  den <- den + (card[i] * card[j])
                }
            }
        }
    }
    distname <- deparse(substitute(dist))
    out <- list(ptc = simptc, ctc = simctc, ratio = (sumnum/num)/(sumden/den), 
                clustering = c, distname=distname, names=attr(dist,'Labels'))
    attr(out,"call") <- call
    attr(out,"class") <- "partana"
    invisible(out)
}

partana.partition <- function (c,dist=NULL)
{
    if (!is.null(c$dist)) {
        tmp <- c$dist
    } else if (!is.null(dist)) {
        tmp <- dist
    } else {
        stop('Your partition object did not contain the dissimilarity object, and you did not provide one')
    }
    attr(tmp, "class") <- "dist"
    out <- partana(c$clustering, tmp)
    out
}


partana.stride <- function(c,dist)
{
    if (class(c) != 'stride')
        stop("The first argument must be of class 'stride'")
    res <- rep(NA,ncol(c$clustering))
    for (i in 1:ncol(c$clustering)) {
        res[i] <- partana(c$clustering[,i],dist)$ratio
    }
    clusters <- c$seq
    ratio <- res
    out <- data.frame(clusters,ratio)
    out
}

plot.partana <- function(x,panel='all',zlim=range(x$ptc),col=heat.colors(12), ...)
{
    numclu <- ncol(x$ptc)
    numplt <- nrow(x$ptc)
    set <- matrix(nrow=0,ncol=numclu)
    card <- rep(0,numclu)
    for (i in 1:numclu) {
        card[i] <- sum(x$clustering==i)
    }
    for (i in 1:numclu) {
        if (card[i] > 0) {
            tmp <- x$ptc[x$clustering==i,]
            if (card[i] > 1) {
                tmp <- tmp[rev(order(tmp[,i])),]
            }
            set <- rbind(set,tmp)
        } 
    }
    if (panel == 'all' || panel == 1) {
        image(seq(1:numplt),seq(1:numclu),set,zlim=zlim,col=col,
            main="Plot-to-Set Similarity",xlab="Plots",ylab="Set")      
        if (panel == 'all') 
            readline("Hit return to continue\n")
    }
    if (panel == 'all' || panel == 2) {
        image(seq(1:numclu),seq(1:numclu),x$ctc,zlim=zlim,col=col,
            main="Set-to-Set Similarity",xlab="Set",ylab="Set")
        if ((panel == 'all' || panel ==3) && length(x$ratio) > 1) 
            readline("Hit return to continue\n")
    }
    if ((panel == 'all' || panel ==3) && length(x$ratio) > 1) {
        plot(x$ratio,type='b')
    }
}
 
summary.partana <- function (object, ...) 
{
    cat(paste("Number of clusters = ", nrow(object$ctc), "\n"))
    print(table(object$clustering))
    cat("\n")
    if (nrow(object$ctc) < 11) {
        print(object$ctc)
    }
    else {
        cat("Mean Within-cluster similarities\n\n")
        for (i in 1:nrow(object$ctc)) {
            cat(paste(i, format(object$ctc[i, i], digits = 4), 
                "\n"))
        }
    }
    if (length(object$ratio) > 1) {
        cat(paste("\nRatio of Within-cluster similarity/Among-cluster similarity = ",
             format(object$ratio[object$numitr],digits=4),"in",object$numitr,"iterations\n"))
    }
    else {
        cat(paste("\nRatio of Within-cluster similarity/Among-cluster similarity = ",
            format(object$ratio,digits=4),"\n"))
    }
}

