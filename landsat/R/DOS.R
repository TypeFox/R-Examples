DOS <-
function(sat=5, scattering.coef=c(-4, -2, -1, -.7, -.5), SHV, SHV.band, gain, offset, Grescale, Brescale, sunelev, edist, Esun=c(198.3, 179.6, 153.6, 103.1, 22, 8.34), blackadjust = 0.01)
{

### Improved Dark Object Subtraction method from Chavez 1989
### Implements the 1% adjustment (can be changed with blackadjust argument)

### Dark Object Subtraction method from Chavez 1988
### with some modifications 
### calculates DN to subtract for each band for a range of
### scattering.coef: scattering coefficients. Default values:
###	-4.0: Very Clear	SHV <= 55
###	-2.0: Clear		SHV 56-75
### 	-1.0: Moderate		SHV 76-95
###	-0.7: Hazy		SHV 96-115
###	-0.5: Very Hazy		SHV >115
### sat: 5 or 7 for Landsat platform
### sunelev: sun elevation in degrees
### satzenith: satellite zenith angle in degrees (= 0 for Landsat)
### edist: Earth-Sun distance (calculated from DOY)
### Esun: extrasolar radiation
### gain and offset or gain and bias
### SHV: lowest DN value; from a black object
### SHV.band: band from which SHV was taken
### blackadjust: lowest DN not from an absolutely black object; reduce it by 1% or other value

### returns the mean version of Chavez 1988 and a slightly more
### sophisticated approximation over the entire band's wavelengths

    if(sat == 5)
        bands <- data.frame(
            lmin=c(0.45, 0.52, 0.63, 0.76, 1.55, 2.08),
            lmax=c(0.52, 0.60, 0.69, 0.90, 1.75, 2.35))
    else if(sat == 7)
        bands <- data.frame(
            lmin=c(0.45, 0.52, 0.63, 0.77, 1.55, 2.09),
            lmax=c(0.52, 0.60, 0.69, 0.90, 1.75, 2.35))
    else stop("Unknown satellite.\n")

    rownames(bands) <- c("band1", "band2", "band3", "band4", "band5", "band7")


    ### Chavez 1988 Table 1
    ### Values of specific functions for the Landsat bands

    scattering.mean <- matrix(apply(bands, 1, mean), byrow=FALSE, nrow=nrow(bands), ncol=length(scattering.coef))
    rownames(scattering.mean) <- rownames(bands)
    colnames(scattering.mean) <- paste("coef", scattering.coef, sep="")
    scattering.mean <- sweep(scattering.mean, 2, scattering.coef, "^")
    scattering.mean.pct <- sweep(scattering.mean, 2, apply(scattering.mean, 2, sum), "/")

    # alternate version using curve approximation
    scattering.approx <- matrix(NA, nrow=nrow(bands), ncol=length(scattering.coef))
    rownames(scattering.approx) <- rownames(bands)
    colnames(scattering.approx) <- paste("coef", scattering.coef, sep="")

    grain <- 0.0001
    for(i in 1:nrow(bands)) {
        thisband <- seq(bands[i, 1], bands[i, 2], by=grain)
        for(j in 1:length(scattering.coef)) {
            scattering.approx[i, j] <- mean(thisband ^ scattering.coef[j])
        }
    }
    scattering.approx.pct <- sweep(scattering.approx, 2, apply(scattering.approx, 2, sum), "/")

    ### Chavez 1988 Table 2
    ### Multiplication factors to predict haze values in other spectral 
    ### bands given a starting haze value and band

    corrband.mean <- scattering.mean[SHV.band, ]
    corrband.mean <- sweep(scattering.mean, 2, corrband.mean, "/")

    corrband.approx <- scattering.approx[SHV.band, ]
    corrband.approx <- sweep(scattering.approx, 2, corrband.approx, "/")

    ### Chavez 1988 Table 3; Chavez 1989 Table 1
    ### Gain, offset and normalization factors

    # most new references provide gain and bias
    # need gain and offset

    # most new references provide gain and bias
    # want gain and offset
    if(missing(offset)) {
        offset <- -1 * Brescale / Grescale
        gain <- 1/Grescale
    }

    NORM <- gain / gain[SHV.band]

    ### convert sunelev in degrees to sun zenith angle in radians
    suntheta <- (90-sunelev) * pi / 180
    suntheta <- cos(suntheta)

    ### Calculate Eo from Esun and edist
    Eo <- Esun[SHV.band]/edist^2
        
    # subtract 1% - assume that black in the image is not really black
    SHV <- SHV - gain[SHV.band] * blackadjust * Eo * suntheta / pi


    ## from here on, follow DOS algorithm of Chavez 1988


    SHV <- SHV - offset[SHV.band]

    DNfinal.mean <- SHV * corrband.mean 
    DNfinal.mean <- sweep(DNfinal.mean, 1, NORM, "*")
    DNfinal.mean <- sweep(DNfinal.mean, 1, offset, "+")

    DNfinal.approx <- SHV * corrband.approx 
    DNfinal.approx <- sweep(DNfinal.approx, 1, NORM, "*")
    DNfinal.approx <- sweep(DNfinal.approx, 1, offset, "+")


    list(DNfinal.mean = DNfinal.mean, DNfinal.approx = DNfinal.approx)

}

