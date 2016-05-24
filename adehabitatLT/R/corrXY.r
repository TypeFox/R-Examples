.convtime <- function(step, units)
{
        if (units=="min")
            step <- step*60
        if (units=="hour")
            step <- step*60*60
        if (units=="day")
            step <- step*60*60*24
        return(step)
}


.corrXY <- function(x, y, dab, daa)
{
    x1 <- x[1]
    y1 <- y[1]

    if (daa[1]!=dab[1]) {
        dt <- daa[2]-daa[1]
        dif <- abs(daa[1]-dab[1])
        alpha <- atan2( (y[2]-y[1]), (x[2]-x[1]) )
        r <- sqrt( (x[2]-x[1])^2 + (y[2]-y[1])^2 )*dif/dt

        if (daa[1] > dab[1]) {
            x1 <- x[1]+cos(alpha)*r
            y1 <- y[1]+sin(alpha)*r
        }
        if (daa[1] < dab[1]) {
            alpha <- pi + alpha
            x1 <- x[1]+cos(alpha)*r
            y1 <- y[1]+sin(alpha)*r
        }
    }

    xn <- x[length(y)]
    yn <- y[length(y)]
    if (daa[length(daa)]!=dab[length(dab)]) {

        dt <- daa[length(daa)]-daa[length(daa)-1]
        dif <- abs(daa[length(daa)]-dab[length(daa)])
        alpha <- atan2((y[length(daa)]-y[length(daa)-1]),
                       (x[length(daa)]-x[length(daa)-1]))

        r <- sqrt((x[length(daa)]-x[length(daa)-1])^2 +
                  (y[length(daa)]-y[length(daa)-1])^2)*dif/dt

        if (daa[length(dab)]<dab[length(dab)]) {
            alpha <- pi + alpha
            xn <- x[length(daa)]+cos(alpha)*r
            yn <- y[length(daa)]+sin(alpha)*r
        } else {
            xn <- x[length(daa)]+cos(alpha)*r
            yn <- y[length(daa)]+sin(alpha)*r
        }
    }

    x[1] <- x1
    y[1] <- y1
    x[length(x)] <- xn
    y[length(x)] <- yn

    xb <- x[1:(length(x)-2)]
    xc <- x[2:(length(x)-1)]
    xa <- x[3:length(x)]
    yb <- y[1:(length(x)-2)]
    yc <- y[2:(length(x)-1)]
    ya <- y[3:length(x)]

    daat <- daa
    dabt <- dab

    if (any(daa[-c(1,length(daa))]<dab[-c(1,length(daa))])) {
        dt <- diff(daa[1:(length(daa)-1)])
        daa <- daa[2:(length(x)-1)]
        dab <- dab[2:(length(x)-1)]
        dif <- abs(daa-dab)

        alpha <- atan2( (yc-yb),(xc-xb) )
        r <- sqrt( (xc-xb)^2 + (yc-yb)^2 )*dif/dt

        xc[daa<dab] <- xc[daa<dab] + cos(alpha[daa<dab] + pi)*r[daa<dab]
        yc[daa<dab] <- yc[daa<dab] + sin(alpha[daa<dab] + pi)*r[daa<dab]
    }

    daa <- daat
    dab <- dabt

    if (any(daa[-c(1,length(daa))] > dab[-c(1,length(daa))])) {
        dt <- diff(daat[2:(length(daat))])
        daa <- daa[2:(length(x)-1)]
        dab <- dab[2:(length(x)-1)]
        dif <- abs(daa-dab)

        alpha <- atan2( (ya-yc),(xa-xc) )
        r <- sqrt( (xa-xc)^2 + (ya-yc)^2 )*dif/dt
        xc[daa>dab] <- xc[daa>dab] + cos(alpha[daa>dab])*r[daa>dab]
        yc[daa>dab] <- yc[daa>dab] + sin(alpha[daa>dab])*r[daa>dab]
    }

    xn <- c(x1, xc, xn)
    yn <- c(y1, yc, yn)

    return(list(x=xn,y=yn))
}

