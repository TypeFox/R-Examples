############################################################################################
## package 'secr'
## make.tri.R
## last changed 2010 07 29
############################################################################################

## make.hex <- function (r = 0.5, xmin = 0, xmax = 1, ymin = 0, ymax = 1) {
##     dx1 <- r * sin (pi / 6)
##     dx2 <- r
##     dy1 <- 2 * r * cos (pi / 6)
##     dy2 <- r * cos (pi / 6)
##     nx <- (xmax-xmin) / (2*dx1+dx2) + 2
##     xvals <- cumsum (c (0, rep(c(dx1, dx2), nx)))
##     xvals <- xvals[xvals<=xmax]
##     temp <- expand.grid(x = xvals, y = seq(ymin, ymax, dy1))
##     temp$y <- temp$y + rep(c(1,0,0,1) * (dy1 - dy2), length=length(xvals))
##     temp[temp$y < ymax,]
## }
## eqscplot(make.hex(r=0.1))

hexagon <- function (r = 20, centre = c(0,0)) {
    ## vertices of a closed hexagon
    theta <- seq(-4*pi/6, 9*pi/6, pi/3)
    temp <- data.frame(x= r * cos(theta), y = r * sin(theta))
    sweep(temp, MARGIN = 2, STATS = centre, '+')
}

make.tri <- function (nx=10, ny=12, spacing = 20, detector='multi', originxy=c(0,0))
{
    if (!( detector %in% .localstuff$validdetectors ))
        stop ("invalid detector type")
    dy2 <- spacing * cos (pi / 6)
    xvals <- seq(from = originxy[1], length = nx, by = spacing)
    yvals <- seq(from = originxy[2], length = ny, by = dy2)
    grid <- expand.grid(x = xvals, y = yvals)
    oddrow <- as.logical(rep ((0:(ny-1)) %% 2, rep(nx, ny)))
    grid$x[oddrow] <- grid$x[oddrow] + spacing/2
    attr(grid, 'detector')   <- detector
    attr(grid, 'class')      <- c('traps', 'data.frame')
    attr(grid, 'spacex')     <- NULL
    attr(grid, 'spacey')     <- NULL
    attr(grid, 'spacing')    <- spacing
    attr(grid, 'usage')      <- NULL
    attr(grid, 'covariates') <- NULL
    grid
}

make.hex <- function (nx=10, ny=12, spacing = 20, detector='multi',
    originxy=c(0,0)) {
    dy2 <- spacing * cos (pi / 6)
    xvals <- seq(from = originxy[1], length = nx, by = spacing)
    yvals <- seq(from = originxy[2], length = ny, by = dy2)
    grid <- expand.grid(x = xvals, y = yvals)
    oddrow <- as.logical(rep ((0:(ny-1)) %% 2, rep(nx, ny)))
    third  <- as.logical(rep ((0:(nx-1)) %% 3,ny))
    grid$x[oddrow] <- grid$x[oddrow] + 3*spacing/2
    grid <- grid[third,]
    attr(grid, 'detector')   <- detector
    attr(grid, 'class')      <- c('traps', 'data.frame')
    attr(grid, 'spacex')     <- NULL
    attr(grid, 'spacey')     <- NULL
    attr(grid, 'spacing')    <- spacing
    attr(grid, 'usage')      <- NULL
    attr(grid, 'covariates') <- NULL
    grid
}

ringID <- function (hextraps) {
    centrexy <- apply(hextraps, 2, mean)
    dist <- function(xy) sqrt(sum((xy-centrexy)^2))
    radians <- function(xy) { dxy <- xy-centrexy; atan2(dxy[2], dxy[1]) }
    d <- apply(hextraps, 1, dist)
    angle <- apply(hextraps, 1, radians)
    ring <- trunc(d / attr(hextraps, 'spacing') + 0.5) + 1
    if (max(ring)>26)
        stop ("limited to 26 rings")
    ringorder <- order(ring, angle)
    hextraps[,] <- hextraps[ringorder,]
    ring <- ring[ringorder]
    ## sequence number within ring
    seqnum <- lapply(split(1:nrow(hextraps), ring), function(x) 1:length(x))
    row.names(hextraps) <- paste(LETTERS[ring], unlist(seqnum), sep='')
    hextraps
}

clip.hex <- function (traps, side = 20, centre = c(50, 60*cos(pi/6)),
    fuzz = 1e-3, ID = 'num', ...) {
    hex <- hexagon(side+fuzz, centre)
    hex <- matrix(unlist(hex), ncol = 2)
    OK <- pointsInPolygon(traps, hex)
    temp <- subset(traps, OK, ...)
    if (!is.null(ID)) {
        if (!(ID %in% c('num', 'alpha')))
            stop ("ID must be 'num' or 'alpha'")
        if (ID == 'num') row.names(temp) <- 1:nrow(temp)
        if (ID == 'alpha') temp <- ringID(temp)
    }
    temp
}

## library(MASS)
## library(secr)
## temp <- make.tri(spacing=10, nx=10, ny=12)
## temp <- clip.hex(temp, side=30, ID='num', )
## plot(temp, bor=5,off=1)
## text(temp$x, temp$y, row.names(temp))
##
## temp <- clip.hex(temp, ID='alpha')
## plot(temp, bor=5,off=1)
## text(temp$x, temp$y, row.names(temp))
##
## plot(make.tri(spacing=10))
## points(clip.hex(make.tri(spacing=10), hex.s=30), pch=16)
## points(clip.hex(make.tri(spacing=10), hex.s=50), pch=16)
