topocorr <-
function(x, slope, aspect, sunelev, sunazimuth, method="cosine", na.value=NA, GRASS.aspect=FALSE, IL.epsilon=0.000001)
{
# topographic correction for image x based on
# topography and sun location

# IL.epsilon: if IL == 0, the corrected value is Inf (division by zero)
# adding a tiny increment eliminates the Inf

## aspect may be GRASS output: counterclockwise from east
## or nonGRASS output: clockwise from north
## require the latter for further calculations
## because sunazimuth is always measured clockwise from north
    if(GRASS.aspect) {
        aspect <- as.matrix(aspect)
        aspect <- -1 * aspect + 90
        aspect <- (aspect + 360) %% 360
    }

# all inputs are in degrees, but we need radians
    slope <- (pi/180) * as.matrix(slope)
    aspect <- (pi/180) * as.matrix(aspect)
    sunzenith <- (pi/180) * (90 - sunelev)
    sunazimuth <- (pi/180) * sunazimuth

    x.orig <- x
    x <- as.matrix(x)
    x[x == na.value] <- NA

    IL <- cos(slope) * cos(sunzenith) + sin(slope) * sin(sunzenith) * cos(sunazimuth - aspect)
    IL[IL == 0] <- IL.epsilon

        METHODS <- c("cosine", "improvedcosine", "minnaert", "minslope", "ccorrection", "gamma", "SCS", "illumination")
        method <- pmatch(method, METHODS)
        if (is.na(method)) 
            stop("invalid method")
        if (method == -1) 
            stop("ambiguous method")

    if(method == 1){
        ## Cosine method
        xout <- x * (cos(sunzenith)/IL)
    }
    else if(method == 2) {
    ## Improved cosine method
        ILmean <- mean(as.vector(IL), na.rm=TRUE)
        xout <- x + (x * (ILmean - IL)/ILmean)
    }
    else if(method == 3) {
        ## Minnaert
        ## K is between 0 and 1
        ## only use points with greater than 5% slope
        targetslope <- atan(.05)

        if(all(x[slope >= targetslope] < 0, na.rm=TRUE)) {
            K <- 1
        }
        else {
            # IL can be <=0 under certain conditions
            # but that makes it impossible to take log10 so remove those elements
            K <- data.frame(y = as.vector(x[slope >= targetslope]), x = as.vector(IL[slope >= targetslope])/cos(sunzenith))
            K <- K[!apply(K, 1, function(x)any(is.na(x))),]
            K <- K[K$x > 0, ]
            K <- K[K$y > 0, ]

            K <- lm(log10(K$y) ~ log10(K$x))
            K <- coefficients(K)[[2]] # need slope
            if(K > 1) K <- 1
            if(K < 0) K <- 0
        }

        xout <- x * (cos(sunzenith)/IL) ^ K
    }
    else if(method == 4) {
        ## Minnaert with slope
        ## K is between 0 and 1
        ## only use points with greater than 5% slope
        targetslope <- atan(.05)

        if(all(x[slope >= targetslope] < 0, na.rm=TRUE)) {
            K <- 1
        }
        else {
            # IL can be <=0 under certain conditions
            # but that makes it impossible to take log10 so remove those elements
            K <- data.frame(y=as.vector(x[slope >= targetslope]), x=as.vector(IL[slope >= targetslope])/cos(sunzenith))
            K <- K[!apply(K, 1, function(x)any(is.na(x))),]
            K <- K[K$x > 0, ]
            K <- K[K$y > 0, ]

            K <- lm(log10(K$y) ~ log10(K$x))
            K <- coefficients(K)[[2]] # need slope
            if(K > 1) K <- 1
            if(K < 0) K <- 0
        }

        xout <- x * cos(slope) * (cos(sunzenith) / (IL * cos(slope))) ^ K
    }
    else if(method == 5) {
        ## C correction
        band.lm <- lm(as.vector(x) ~ as.vector(IL))
        C <- coefficients(band.lm)[[1]]/coefficients(band.lm)[[2]]

        xout <- x * (cos(sunzenith) + C) / (IL + C)
    }
    else if(method == 6) {
        ## Gamma
        ## assumes zenith viewing angle
        viewterrain <- pi/2 - slope
        xout <- x * (cos(sunzenith) + cos(pi / 2)) / (IL + cos(viewterrain))
    }
    else if(method == 7) {
        ## SCS method from GZ2009
        xout <- x * (cos(sunzenith) * cos(slope))/IL
    }
    else if(method == 8) {
        ## illumination only
        xout <- IL
    }

    ## if slope is zero, reflectance does not change
    if(method != 8) 
        xout[slope == 0 & !is.na(slope)] <- x[slope == 0 & !is.na(slope)]

    ## if x was a SpatialGridDataFrame, return an object of the same class
    if(class(x.orig) == "SpatialGridDataFrame") {
        x.orig@data[,1] <- as.vector(xout)
        xout <- x.orig
    }


    xout

}

