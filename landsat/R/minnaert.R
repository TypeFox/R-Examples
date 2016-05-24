
minnaert <- function(x, slope, aspect, sunelev, sunazimuth, na.value=NA, GRASS.aspect=FALSE, IL.epsilon=0.000001, slopeclass = c(1, 5, 10, 15, 20, 25, 30, 45), coverclass)
{
# pixel-based Minnaert correction
# Lu et al. 2008. PE&RS 74:1343-1350.
# topographic correction for image x based on
# topography and sun location
# require(mgcv)

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
    sloper <- (pi/180) * as.matrix(slope)
    sloped <- as.matrix(slope)
    aspect <- (pi/180) * as.matrix(aspect)
    sunzenith <- (pi/180) * (90 - sunelev)
    sunazimuth <- (pi/180) * sunazimuth

    x.orig <- x
    x <- as.matrix(x)
    x[x == na.value] <- NA

    IL <- cos(sloper) * cos(sunzenith) + sin(sloper) * sin(sunzenith) * cos(sunazimuth - aspect)
    IL[IL == 0] <- IL.epsilon

    if(missing(coverclass)) 
        coverclass <- rep(TRUE, length(as.vector(x)))

    ## Minnaert
    ## K is between 0 and 1
        # IL can be <=0 under certain conditions
        # but that makes it impossible to take log10 so remove those elements
    K <- data.frame(x = as.vector(x), IL = as.vector(IL), slope=as.vector(sloped))
    K <- K[coverclass, ]
    K <- K[!apply(K, 1, function(x)any(is.na(x))),]
    K <- K[K$x > 0, ]
    K <- K[K$IL > 0, ]

    ## calculate overall K value; only use points with greater than 5% slope
    targetslope <- (180/pi) * atan(.05)
    allcoef <- coefficients(lm(log10(K$x)[K$slope >= targetslope] ~ log10(K$IL/cos(sunzenith))[K$slope >= targetslope]))[[2]]

    results <- data.frame(matrix(0, nrow=length(slopeclass)-1, ncol=3))
    colnames(results) <- c("midpoint", "n", "k")
    results[,1] <- diff(slopeclass)/2 + slopeclass[1:length(slopeclass)-1]

    K.cut <- as.numeric(cut(K$slope, slopeclass)) # don't use slopes outside slopeclass range
    if(nrow(results) != length(table(K.cut))) stop("slopeclass is inappropriate for these data (empty classes)\n")
    results[,2] <- table(K.cut)

    #            sapply(unique(K.cut), function(i)coefficients(lm(log10(K$x)[K.cut == i] ~ log10(K$IL/cos(sunzenith))[K.cut == i]))[[2]])
   for(i in sort(unique(K.cut[!is.na(K.cut)]))) {
        results[i, 3] <- coefficients(lm(log10(K$x)[K.cut == i] ~ log10(K$IL/cos(sunzenith))[K.cut == i]))[[2]]
    }

    model <- with(results, gam(k ~ s(midpoint, k=length(midpoint)-1)))

    K.all <- data.frame(midpoint = as.vector(as.matrix(slope)))
    K.all[K.all > max(slopeclass)] <- max(slopeclass) # if slope is greater than modeled range, use maximum of modeled range
    K.all[K.all < min(slopeclass)] <- 0 # if slope is less than modeled range, treat it as flat
    K.all <- predict(model, newdata=K.all)
    K.all[K.all > 1] <- 1
    K.all[K.all < 0] <- 0

    xout <- as.vector(as.matrix(x)) * (cos(sunzenith)/as.vector(as.matrix(IL))) ^ K.all
    xout[K.all == 0 & !is.na(K.all)] <- as.vector(as.matrix(x))[K.all == 0 & !is.na(K.all)] # don't correct flat areas

    ## if x was a SpatialGridDataFrame, return an object of the same class
    if(class(x.orig) == "SpatialGridDataFrame") {
        x.orig@data[,1] <- as.vector(xout)
        xout <- x.orig
    }

    list(allcoef=allcoef, classcoef=results, model=model, minnaert=xout)

}

