#
# Fit a spherical smoothing spline model
#

smooth.sspline <- function(lon, lat, y, m = 2, smth = 0, lambda = 0)
{
    this.call <- match.call()
    
    if ((m<1)||(m>10)) stop("The order of smoothing should be between 1 and 10!")
    
    if ((smth=1)&&(lambda<0)) stop("The smoothing parameter lambda must be positive!")
    
    n <- length(lon)

    #
    # create the input/ouput/workspace variables.
    #
    varhat <- numeric(1)
    gcv    <- numeric(1)
    coef.d <- numeric(1)
    coef.c <- numeric(n)
    yhat   <- numeric(n)
    S      <- numeric(n)
    tau    <- numeric(n)
    d      <- numeric(n)
    e      <- numeric(n)
    work   <- numeric(8*n)

    #
    # convert the longitudes and latitudes into radiant measure.
    #
    lon <- (lon*pi)/180;  lat <- (lat*pi)/180
    
    #
    # call the Fortran subroutine to calculate the spherical spline.
    #    
    result <- .Fortran("sspline", as.double(lon), as.double(lat),
        as.double(y), as.integer(n), as.integer(m), as.integer(smth),
        as.double(lambda), as.double(gcv), as.double(varhat),
        as.double(coef.c), as.double(coef.d), as.double(yhat),
        as.double(numeric(n^2)), as.integer(n), as.double(S),
        as.double(tau), as.double(work), as.integer(8*n), as.double(d),
        as.double(e), PACKAGE="sspline")[c(7:12)]
    
    #
    # construct a smooth.sspline object from the returned list.
    #
    object <- list(lon = lon, lat = lat, obs = y, m = m,
        smth = smth, lambda = result[[1]], gcv = result[[2]],
        varhat = result[[3]], c = result[[4]], d = result[[5]],
        yhat = result[[6]], call = this.call)

    class(object) <- "smooth.sspline"
    object
}


#
# Define the print method for the class smooth.sspline
#

print.smooth.sspline <- function (x, ...)
{    
    if (!is.null(cl <- x$call))
    {        
        cat("\nCall:\n  ")
        dput(cl)
    }
    cat("\nSample Size n:", length(x$lon), "\n")
    cat("Order of Smooth:", format(x$m), "\n")
    cat("Smoothing Par:", format(x$lambda), "\n")
    cat("GCV Criteria:", format(x$gcv), "\n")
    cat("Estimated Var:", format(x$varhat), "\n\n")
    invisible(x)
}


#
# Define the summary method for the class smooth.sspline
#

summary.smooth.sspline <- function(object, ...)
{
    print(object)
    cat("Coefficient c:\n")
    print(object$c)
    cat("\nCoefficient d:", format(object$d), "\n\n")
}


#
# Define the predict method for the class smooth.sspline
#

predict.smooth.sspline <- function(object, lon, lat, grid = FALSE, ...)
{
    if (missing(lon) && missing(lat))
        return(object$yhat)
    else if (missing(lon))
        lon <- seq(-180, 180, len=101)[2:100]
    else if (missing(lat))
        lat <- seq(-90, 90, len=51)[2:50]
    else   # both lon and lat are NOT missing
        if (!grid && length(lon) != length(lat))
            stop("The longitudes and latitudes must have the same length!\n")
        else if (!grid)
            newdata <- data.frame(lon=lon, lat=lat)
        else
        {
            newdata <- expand.grid(lat, lon)[ , c(2, 1)]
            names(newdata) <- c("lon", "lat")
        }

    newdata$lon <- newdata$lon*pi/180; newdata$lat <- newdata$lat*pi/180
    npred <- nrow(newdata)

    predictor <- .Fortran("ssplfit", as.double(object$c), as.double(object$d),
        as.integer(length(object$c)), as.integer(object$m), 
	as.double(object$lon), as.double(object$lat), as.double(newdata$lon), 
	as.double(newdata$lat), as.integer(npred),
	as.double(numeric(npred)), PACKAGE="sspline")[[10]]
    
    if (grid)
        matrix(predictor, byrow = TRUE, nrow=length(lon))
    else
        predictor
}


#
# define the plot method for the class smooth.sspline
#

plot.smooth.sspline <- function(x, lon, lat, main = "", xlab = "Longitude",
    ylab = "Latitude", key.title = "Temp\n(deg)", ...)
{
    #
    # make the prediction on the given (lon, lat) grid
    #
    predmat <- predict(x, lon = lon, lat = lat, grid = TRUE)
    
    #
    # open a graphics device and plot the filled contour with the world map
    #
    # X11(width = 11, height = 8)
    
    x.at <- seq(-5*pi/6, 5*pi/6, by = pi/6)
    y.at <- seq(-pi/3, pi/3, by = pi/6)
    x.lab <- seq(-150, 150, by = 30)
    y.lab <- seq(-60, 60, by = 30)

    lon <- lon*pi/180; lat <- lat*pi/180
    filled.contour(lon, lat, predmat, color.palette = terrain.colors,
        xlim = c(-pi, pi), ylim = c(-pi/2, pi/2),
        plot.title = title(main = main, xlab = xlab, ylab = ylab),
        plot.axes ={ axis(1, at = x.at, x.lab)
                     axis(2, at = y.at, y.lab)
                     map.world(add = TRUE) },
        key.title = title(main = key.title),
        key.axes  = axis(4, pretty(range(predmat), n=8)),
        lwd = 0.1, ...)

    invisible()    
}

