############################################################################################
## package 'secrlinear'
## make.line.R
## last changed 2014-08-28, 2014-09-04, 2014-09-05
############################################################################################

make.line <- function (SLDF, n = 10, startbuffer = 0, by = 20, endbuffer = 0, cluster = NULL,
                       type = c('fixedstart', 'randomstart', 'centred'), detector = "multi")
{
    alongline <- function(i) {
        len <- lgth[i] - endbuffer
        if (len < startbuffer)
            return (matrix(nrow = 0, ncol = 2))
        if (type %in% c('centred', 'randomstart')) {
            fullline <- (n-1) * by 
            if (!is.null(cluster))
                fullline <- fullline + max(cluster)
            if (type == 'centred') {
                if (fullline < len) {
                    from <- (len-fullline)/2
                }
                else
                    from <- by/2
            }
            else {
                from <- (len-fullline) * runif(1)
            }
            len <- min(from + fullline, len)            
        }            
        else from <- startbuffer
        along <- seq(from = from, by = by, to = len)
        along <- along[along>=0]        
        if (length(along)>n) along <- along[1:n]
        if (!is.null(cluster)) {
            ## replace each point by a cluster
            along <- as.numeric(outer(cluster, along, FUN='+'))
            ## not yet used
            clusterID <- rep(1:length(along), each = length(cluster))
            detectorID <- rep(1:length(cluster), length(along))
        }
        OK <- along <= lgth[i]
        along <- along[OK]        

        xy <- do.call(rbind, xy[[i]])
        dx <- diff(xy[, 1])
        dy <- diff(xy[, 2])
        d <- sqrt(dx^2 + dy^2)
        theta <- atan2(dy, dx)
        cumd <- cumsum(d)
        startvert <- sapply(along, function(x) which.min(cumd <= x))
        if (length(startvert) == 0)
            return (matrix(nrow=0, ncol=2))
        cumd <- c(0, cumd)
        r <- along - cumd[startvert]
        xy1 <- xy[startvert, ]
        px <- r * cos(theta[startvert])
        py <- r * sin(theta[startvert])
        xy1 + cbind(x = px, y = py)
    }
    
    type = match.arg(type)
    if (!is.null(cluster))
        cluster <- cluster - min(cluster)
    if (inherits(SLDF, "linearmask"))
        SLDF <- attr(SLDF, "SLDF")
    lgth <- SpatialLinesLengths(SLDF)
    xy <- coordinates(SLDF)
    tmp <- lapply(1:nrow(SLDF), alongline)
    trps <- do.call(rbind, tmp)
    trps <- data.frame(trps)
    class(trps) <- c("traps", "data.frame")
    detector(trps) <- detector
    spacing(trps) <- by
    attr(trps, "spacex") <- by
    attr(trps, "spacey") <- 0
    attr(trps, "SLDF") <- SLDF
    trps
}

## alongmask <- function (mask, pos)
## {
    ## SLDF <- attr(mask, "SLDF")
    ## lgth <- SpatialLinesLengths(SLDF)
    ## xy <- coordinates(SLDF)
    ## tmp <- lapply(1:nrow(SLDF), alongline, along = pos, lgth, xy)
    ## xy2 <- do.call(rbind, tmp)
    ## data.frame(xy2)
## }