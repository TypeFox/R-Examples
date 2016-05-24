
".speedStats" <- function(x, vdist)
{
    ## Value: A 3-column matrix with total distance, mean speed, and angle
    ## for a section of a dive
    ## --------------------------------------------------------------------
    ## Arguments: x=matrix with a dive's section data (time, speed);
    ## vdist=vertical distance travelled during this time If vdist is
    ## missing, then it's all horizontal movements (no angles)
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    if (nrow(x) > 1) {
        speed <- x[-1, 2]
        time <- x[, 1]
        difft <- diff(as.numeric(time)) # seconds
        mspeed <- mean(speed, na.rm=TRUE)
        tdist <- sum(difft * speed, na.rm=TRUE)
        if (!missing(vdist)) {
            angle <- asin(ifelse(vdist < tdist, vdist/tdist, NA)) * (180 / pi)
            cbind(tdist=tdist, mean.speed=mspeed, angle=angle)
        } else {
            cbind(tdist=tdist, mean.speed=mspeed, angle=NA)
        }
    } else {
        matrix(ncol=3)
    }
}
