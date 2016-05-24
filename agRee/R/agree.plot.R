panel.identity <- function(x, y, ...){
    rr <- range(c(x,y))
    rr[1] <- rr[1] - 0.10*diff(rr)
    rr[2] <- rr[2] + 0.10*diff(rr)
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(rr[1], rr[2], rr[1], rr[2]))
    points(x, y)
    abline(a=0, b=1, lwd=2)
}

panel.blandaltman <- function(x, y, ...){
    mean <- (y+x) / 2
    diff <- y - x
    mu.diff <- mean(diff)
    sd.diff <- sd(diff)
    l.limit <- mu.diff - 3*sd.diff
    u.limit <- mu.diff + 3*sd.diff
    
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(min(min(x), range(mean)[1]),
                max(max(x), range(mean)[2]),
                min(min(y)-0.10*diff(range(y)), l.limit),
                max(max(y)+0.10*diff(range(y)), u.limit)))
    
    lines(mean, diff, type="p")
    abline(h=0, lwd=2)
    abline(h=mu.diff, col="red", lty=3, lwd=2)
    abline(h=mu.diff + 2*sd.diff, lty=3, lwd=2, col="red")
    abline(h=mu.diff - 2*sd.diff, lty=3, lwd=2, col="red")
}

agree.plot <- function(ratings, NAaction=c("fail", "omit")){
    if(!is.matrix(ratings) || ncol(ratings) < 2 || nrow(ratings) < 2)
        stop("'ratings' has to be a matrix of at least two columns and two rows.")
    na <- match.arg(NAaction)
    ratings <- switch(na,
                      fail = na.fail(ratings),
                      omit = na.omit(ratings))
    if(!is.matrix(ratings) || ncol(ratings) < 2|| nrow(ratings) < 2)
        stop("'ratings' has to be a matrix of at least two columns and two rows after removing missing values.")

    pairs(ratings, lower.panel=panel.blandaltman, upper.panel=panel.identity,
          xaxt="n", yaxt="n", labels=paste("Rater", 1:ncol(ratings)))

}



