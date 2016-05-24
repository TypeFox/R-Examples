radiocorr <-
function(x, gain, offset, Grescale, Brescale, sunelev, satzenith=0, edist, Esun, Lhaze, method="apparentreflectance")
{

### Radiometric correction following one of four models:
###   - 1: Absolute Radiance
###   - 2: Dark Object Subtraction (DOS) of Chavez 1989
###   - 3: COSTZ of Chavez 1996
###   - 4: DOS4 of SWS+2001
### x: image data
### gain and offset or gain and bias or Lmin and Lmax
### sunelev: sun elevation in degrees
### satzenith: satellite zenith angle in degrees (= 0 for Landsat)
### edist: Earth-Sun distance (calculated from DOY)
### Esun: extrasolar radiation
### Lhaze: result of DOS() for this band, if needed
### crosscalib: cross-calibration coefficient of Teillet et al. 2006 or other
###    multiplicative adjustment to be applied to at-sensor reflectance - NOT DONE!

    results <- x
    x <- as.vector(as.matrix(x))

    METHODS <- c("apparentreflectance", "DOS", "COSTZ", "DOS4")
    method <- pmatch(method, METHODS)
    if (is.na(method)) 
        stop("invalid method")
    if (method == -1) 
        stop("ambiguous method")

    suntheta <- (90-sunelev) * pi / 180
    suntheta <- cos(suntheta)

    satzenith <- satzenith * pi / 180
    satphi <- cos(satzenith)

    # most new references provide gain and bias
    # want gain and offset
    if(missing(offset)) {
        offset <- -1 * Brescale / Grescale
        gain <- 1/Grescale
    }

    ### Done with basic setup.

    if(method == 1) {
    ## 1. Apparent Reflectance
        TAUz <- 1.0
        TAUv <- 1.0
        Edown <- 0.0 
        Lhaze <- 0.0
    }
    else if(method == 2) {
    ## 2. DOS
        TAUz <- 1.0
        TAUv <- 1.0
        Edown <- 0.0 
        if(missing(Lhaze)) stop("This model requires Lhaze to be specified.\n")	
    }
    else if(method == 3) {
    ## 3. COSTZ
        TAUz <- suntheta
        TAUv <- satphi
        Edown <- 0.0
        if(missing(Lhaze)) stop("This model requires Lhaze to be specified.\n")
    } 
    else if(method == 4) {
    ## 4. DOS4 of SWS+2001
        TAUv <- TAUz <- 1
        taudiff <- 1
            tau <- 9999
            Edown <- 0
            
        Lhaze.orig <- Lhaze
        
        while(abs(taudiff) > 0.0000001) {
            taudiff <- tau
            
            ## if Lhaze is too large, the formula tries to take log of a negative number
            ## iteratively adjust Lhaze downward until it works
            ## This is a lazy kludge!!!

            Eo <- Esun/edist^2

            Lp <- (Lhaze - offset) / gain - 0.01 * (Eo * suntheta * TAUz + Edown) * TAUv / pi

            taustep <- 1 - (4 * pi * Lp) / (Eo * suntheta)

            while(taustep < 0) {
                Lhaze <- Lhaze - 1
                Lp <- (Lhaze - offset) / gain - 0.01 * (Eo * suntheta * TAUz + Edown) * TAUv / pi
                taustep <- 1 - (4 * pi * Lp) / (Eo * suntheta)
            }
            
            tau <- -1 * suntheta * log(1 - (4 * pi * Lp) / (Eo * suntheta))
            TAUv <- exp(-1 * tau / satphi)
            TAUz <- exp(-1 * tau / suntheta)		
            Edown <- pi * Lp
            
                    taudiff <- taudiff - tau

        }
            
        if(!identical(Lhaze.orig, Lhaze)) warning(paste("Lhaze adjusted from ", Lhaze.orig, " to ", Lhaze, sep=""))
        
        if(missing(Lhaze)) stop("This model requires Lhaze to be specified.\n")
    #-#

    }


    ## Applying the models
    ## REF <-  (pi * edist^2 * (Lsat - Lhaze)) / (TAUv * (Esun * cos(sunzenith) * TAUz + Edown))

    ## First convert DN to at-sensor radiance
    ## Lhaze output from DOS() is in DN, so this is done as a separate step

    # subtract Lhaze and convert DN to radiance
    x <- x - Lhaze
    x <- (x - offset) / gain

    ## proceed with radiometric correction
    ## calculate at-surface reflectance

    x <-  (pi * edist^2 * x) / (TAUv * (Esun * suntheta * TAUz + Edown))

    # return the same structure as the input values
    if(class(results) == "SpatialGridDataFrame")
        results@data[,1] <- x
    else if(is.data.frame(x))
        results <- data.frame(matrix(x, nrow=nrow(results), ncol=ncol(results)))
    else # return a matrix 
        results <- x
    
    results
}

