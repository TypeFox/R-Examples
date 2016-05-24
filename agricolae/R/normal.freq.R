`normal.freq` <-
function (histogram,frequency=1, ...)
{
xx <- histogram$mids
if(frequency==1)yy<-histogram$counts
if(frequency==2)yy<-histogram$counts/sum(histogram$counts)
if(frequency==3)yy<-histogram$density   
    media <- sum(yy * xx)/sum(yy)
    variancia <- sum(yy * (xx - media)^2)/sum(yy)
    zz <- histogram$breaks
    x1 <- xx[1] - 4 * (zz[2] - zz[1])
    z <- length(zz)
    x2 <- xx[z - 1] + 4 * (zz[z] - zz[z - 1])
    x <- seq(x1, x2, length = 200)
    y <- rep(0, 200)
    area <- 0
    for (k in 1:(z - 1)) area = area + yy[k] * (zz[k + 1] - zz[k])
    for (i in 1:200) {
        y[i] <- area * exp(-((x[i] - media)^2)/(2 * variancia))/sqrt(2 * 
            pi * variancia)
    }
    lines(x, y, ...)
    abline(h = 0)
}

