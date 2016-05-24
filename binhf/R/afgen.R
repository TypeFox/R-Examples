"afgen" <-
function (xgrid = seq(0, 1, length = 21), ygrid = seq(0, 1, length = 21), 
    samples = 1000, binsize = 32) 
{
    af1vals <- array(0, c(length(xgrid), length(ygrid), samples))
    ansc1vals <- array(0, c(length(xgrid), length(ygrid), samples))
    free1vals<-ansc1vals
    for (i in 1:length(xgrid)) {
        for (j in 1:length(ygrid)) {
            x1vals <- rbinom(samples, binsize, xgrid[i])
            x2vals <- rbinom(samples, binsize, ygrid[j])
            x1x2 <- rbind(x1vals, x2vals)
            maxi <- which.max(c(xgrid[i], ygrid[j]))
            ansc1vals[i, j, ] <- ansc(x1x2[maxi, ], binsize)
		free1vals[i, j, ] <- free(x1x2[maxi, ], binsize)
            af1vals[i, j, ] <- (x2vals - x1vals)/sqrt((x1vals + 
                x2vals) * ((2 * binsize) - (x1vals + x2vals))/(2 * 
                binsize))
        }
    }
    af1vals[which(abs(af1vals) == Inf)] <- 0
    af1vals[which(is.na(af1vals))] <- 0
    ansc1vals[which(abs(ansc1vals) == Inf)] <- 0
    ansc1vals[which(is.na(ansc1vals))] <- 0
    free1vals[which(abs(free1vals) == Inf)] <- 0
    free1vals[which(is.na(free1vals))] <- 0

    return(list(a = af1vals, b = ansc1vals, c = free1vals))
}

