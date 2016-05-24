tasscap <-
function(basename, sat=7)

{
    # basename is the name of the band data files, which will have the band number appended
    # should be in at-sensor reflectance (rc1)
    # sat: 5 = Landsat 5 (TM) or 7 = Landsat 7 (ETM+)
        
        # original papers
        # Kauth and Thomas
        # Crist and Cicone

    band1 <- get(paste(basename, "1", sep=""))
    band2 <- get(paste(basename, "2", sep=""))
    band3 <- get(paste(basename, "3", sep=""))
    band4 <- get(paste(basename, "4", sep=""))
    band5 <- get(paste(basename, "5", sep=""))
    band7 <- get(paste(basename, "7", sep=""))

    if(class(band1) == "SpatialGridDataFrame") {
    	output.sgdf <- band1
	    use.sgdf <- TRUE
    	band1 <- band1@data[,1]
    	band2 <- band2@data[,1]
    	band3 <- band3@data[,1]
    	band4 <- band4@data[,1]
    	band5 <- band5@data[,1]
    	band7 <- band7@data[,1]
    }

    all.bands <- cbind(band1, band2, band3, band4, band5, band7)

        
    if(sat == 7) {
        tc.coef <- matrix(c(
        # Tasseled cap coefficients for Landsat 7 ETM+ at-satellite reflectance from HWY+2002
    # Band 1     Band 2       Band 3     Band 4     Band 5        Band 7     Index
     0.3561,     0.3972,      0.3904,    0.6966,    0.2286,       0.1596,    #  Brightness       
    -0.3344,    -0.3544,     -0.4556,    0.6966,   -0.0242,      -0.2630,    #  Greenness       
     0.2626,     0.2141,      0.0926,    0.0656,   -0.7629,      -0.5388,    #  Wetness          
     0.0805,    -0.0498,      0.1950,   -0.1327,    0.5752,      -0.7775,    #  Fourth           
    -0.7252,    -0.0202,      0.6683,    0.0631,   -0.1494,      -0.0274,    #  Fifth           
     0.4000,    -0.8172,      0.3832,    0.0602,   -0.1095,       0.0985     #  Sixth            
    ), ncol=6, byrow=TRUE)
    } else if(sat == 5) {
        tc.coef <- matrix(c(
        # TM Tasseled Cap Equivalent Transformation Matrix for Band Reflectance Factor from Crist1985
    # Band 1     Band 2       Band 3     Band 4     Band 5        Band 7     Index
     0.2043,     0.4158,      0.5524,    0.5741,    0.3124,       0.2303,    #  Brightness
    -0.1603,    -0.2819,     -0.4934,    0.7940,    0.0002,      -0.1446,    #  Greenness
     0.0315,     0.2021,      0.3102,    0.1594,    0.6806,      -0.6109,    #  Wetness
    -0.2117,    -0.0284,      0.1302,   -0.1007,    0.6529,      -0.7078,    #  Fourth
    -0.8669,    -0.1835,      0.3856,    0.0408,    0.1132,       0.2272,    #  Fifth
     0.3677,    -0.8200,      0.4354,    0.0518,    0.0066,      -0.0104     #  Sixth
    ), ncol=6, byrow=TRUE)
    } else {
        stop("sat not recognized.\n")
    }

    colnames(tc.coef) <- c("band1", "band2", "band3", "band4", "band5", "band7")
    rownames(tc.coef) <- c("Brightness", "Greenness", "Wetness", "Fourth", "Fifth", "Sixth")
    tc.coef <- t(tc.coef)

    output <- all.bands %*% tc.coef
    output <- as.data.frame(output[,1:3])

    if(use.sgdf) {
    	Brightness <- output.sgdf
	Brightness@data[,1] <- output[, "Brightness"]
    	Greenness <- output.sgdf
	Greenness@data[,1] <- output[, "Greenness"]
    	Wetness <- output.sgdf
	Wetness@data[,1] <- output[, "Wetness"]
	output <- list(Brightness=Brightness, Greenness=Greenness, Wetness=Wetness)
    }

    output
}

