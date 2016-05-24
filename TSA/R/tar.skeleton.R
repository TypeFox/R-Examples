`tar.skeleton` <-
function (object, Phi1, Phi2, thd, d, p,ntransient = 500, n = 500, 
xstart,plot=TRUE,n.skeleton=50) 
{

tar1.skeleton <-
function (x, Phi1, Phi2, thd, d) 
{
    if (x[d] <= thd) {
        newx <- sum(c(1, x) * Phi1)
    }
    else {
        newx <- sum(c(1, x) * Phi2)
    }
    newx
}


if(!missing(object)) {constant1=object$qr1$coefficients[1]
p1=object$p1
p2=object$p2
d=object$d
p=max(p1,p2)
mpd=max(p,d)
Phi1=rep(0,mpd+1)
Phi2=rep(0,mpd+1)
if(object$is.constant1) Phi1a=object$qr1$coefficients else Phi1=c(0,object$qr1$coefficients)
if(object$is.constant2) Phi2a=object$qr2$coefficients else Phi2=c(0,object$qr2$coefficients)
Phi1[seq(Phi1a)]=Phi1a
Phi2[seq(Phi2a)]=Phi2a
thd=object$thd

if(missing(xstart)) xstart=rev(rev(object$y)[1:mpd])
}
if(missing(object) & missing(xstart)) xstart=rep(0,mpd)
        x <- rep(0,ntransient + n)
    x[1:mpd] <- xstart
    for (i in ((mpd + 1):(ntransient + n))) {
        x[i] <- tar1.skeleton(x[(i - 1):(i - mpd)], Phi1 = Phi1, 
            Phi2 = Phi2, thd = thd, d = d)
    }
    tail <- rev(rev(x)[1:n.skeleton])
    nocycle <- TRUE
    i <- 1
    m <- length(tail)
    if (abs(tail[m]) > 10^7) {
        cat("\n skeleton explodes to infinity\n")
        return(invisible())
    }
    while (nocycle & (i < floor(m/2))) {
        nocycle <- any(abs(tail[-(1:i)] - tail[1:(m - i)]) > 
            1e-04)
        if (!nocycle) 
            break
        i <- i + 1
    }
    if (nocycle) {
        cat("\n no limit cycle \n")
        cat(" \n tail part of the skeleton: \n ")
        print(signif(tail, 3))
    }
    else cat(" limit cycle of length ", i, " and it consists of ", 
        signif(limit <- tail[1:i], 5), "\n")
if(plot) plot(y=tail,x=seq(tail),ylab='Skeleton',xlab='t',type='b')
    invisible(tail)
}

