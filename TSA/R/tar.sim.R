`tar.sim` <-
function (object,ntransient = 500, n = 500,  Phi1, 
    Phi2, thd, d, p, sigma1, sigma2, xstart = rep(0, max(p,d)), e) 
{
if(!missing(object)) {

p1=object$p1
p2=object$p2
p=max(p1,p2)
Phi1=rep(0,p+1)
Phi2=rep(0,p+1)
if(object$is.constant1) Phi1[1:(p1+1)]=object$qr1$coefficients  else Phi1[2:(p1+1)]=object$qr1$coefficients 
if(object$is.constant2) Phi2[1:(p2+1)]=object$qr2$coefficients  else Phi2[2:(p2+1)]=object$qr2$coefficients 
d=object$d
thd=object$thd
sigma1=object$rms1^0.5
sigma2=object$rms2^0.5
}
if(missing(xstart)) xstart=rep(0,max(p,d))

    large <- 10^7
    mm=length(xstart)
    x <- rep(0, (mm+ntransient + n))
    boot.start <- length(xstart) + 1
    x[seq(xstart)] <- xstart
    if (missing(e)) 
        e <- rnorm(mm+ntransient + n)
    for (i in (boot.start:(mm+ntransient + n))) {
        if (x[i - d] <= thd) 
            x[i] <-  sum(Phi1 * c(1,x[(i - 1):(i - p)])) + 
                sigma1 * e[i]
        else x[i] <- sum(Phi2 * c(1,x[(i - 1):(i - p)])) + 
            sigma2 * e[i]
    }
    invisible(list(y = x[(ntransient+mm + 1):(ntransient +mm+ n)], 
        e = e[(ntransient +mm+ 1):(ntransient +mm+ n)]))
}

