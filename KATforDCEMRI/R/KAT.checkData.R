##    KATforDCEMRI: a Kinetic Analysis Tool for DCE-MRI
##    Copyright 2014 Genentech, Inc.
##
##    For questions or comments, please contact
##    Gregory Z. Ferl, Ph.D.
##    Genentech, Inc.
##    Development Sciences
##    1 DNA Way, Mail stop 463A
##    South San Francisco, CA, United States of America
##    E-mail: ferl.gregory@gene.com

KAT.checkData <-
    function(file.name, vector.times, map.CC, mask.ROI, vector.AIF, write.data.to.file=TRUE){

        st <- 1
        sp <- dim(mask.ROI)[3]

        time <- as.vector(vector.times)
        cc <- map.CC[,,st:sp,]
        roi <- mask.ROI[,,st:sp]
        aif <- as.vector(vector.AIF)

        ## DETRMINE IF UNITS OF TIME VECTOR ARE IN MINUTES OR SECONDS
        if(max(time)<100)
            t.units <- "minutes"
        if(max(time)>=100)
            t.units <- "seconds"
        if(max(time)>5000 | max(time)<1)
            warning("It looks like the units of your time vector is not minutes or seconds. If this is true then convert to minutes or seconds and try again, otherwise ignore this warning.")
        if(t.units=="minutes")
            time <- time*60

        cat("\n")
        cat("checking dimensions of vectors and arrays...", "\n")
        cat("\n")
        cat("length of vector.times is", length(time), "with units of", t.units, "\n")
        cat("length of vector.AIF is", length(aif), "\n")

        if(length(time) != length(aif))
            stop("vector.times and vector.AIF must be the same length")

        if(length(dim(cc)) != (length(dim(roi))+1))
            stop("map.CC array must have n+1 dimensions when mask.ROI array has n dimensions (one of these arrays may have data corresponding to a single slice while the other has multiple slices)")

        if(length(dim(cc))==4){
            cat("dimensions of map.CC array are", length(cc[,1,1,1]), "x", length(cc[1,,1,1]), "x", length(cc[1,1,,1]), "slices x", length(cc[1,1,1,]), "time points", "\n")

            if(length(time) != length(cc[1,1,1,]))
                stop("length of vector.times and nt dimension of the map.CC array must be equal")

            if(length(aif) != length(cc[1,1,1,]))
                stop("length of vector.AIF and nt dimension of the map.CC array must be equal")

            cat("dimensions of mask.ROI array are", length(roi[,1,1]), "x", length(roi[1,,1]), "x", length(roi[1,1,]), "slices", "\n")

            if(length(roi[,1,1]) != length(cc[,1,1,1]))
                stop("nx dimension of the mask.ROI and map.CC arrays must be equal")

            if(length(roi[1,,1]) != length(cc[1,,1,1]))
                stop("ny dimension of the mask.ROI and map.CC arrays must be equal")

            if(length(roi[1,1,]) != length(cc[1,1,,1]))
                stop("number of slices within mask.ROI and map.CC arrays must be equal")
        }

        if(length(dim(cc))==3){
            cat("dimensions of map.CC array are", length(cc[,1,1]), "x", length(cc[1,,1]), "x", length(cc[1,1,]), "time points", "\n")

            if(length(time) != length(cc[1,1,]))
                stop("length of vector.times and nt dimension of the map.CC array must be equal")

            if(length(aif) != length(cc[1,1,]))
                stop("length of vector.AIF and nt dimension of the map.CC array must be equal")

            cat("dimensions of mask.ROI array are", length(roi[,1]), "x", length(roi[1,]), "\n")

            if(length(roi[,1]) != length(cc[,1,1]))
                stop("nx dimension of the mask.ROI and map.CC arrays must be equal")

            if(length(roi[1,]) != length(cc[1,,1]))
                stop("ny dimension of the mask.ROI and map.CC arrays must be equal")
        }

        cat("\n")
        cat("...vector and array dimensions are okay.", "\n")
        cat("\n")

        if(write.data.to.file==TRUE)
            cat("Saving data in a single R file...", "\n")

        ##file.name <- paste(file.name, "_s", st, "-s", sp, ".RData", sep="")
        file.name <- paste(file.name, ".RData", sep="")

        dcemri.data <- list(time, cc, roi, aif, write.data.to.file)
        names(dcemri.data) <- c("vectorTimes", "mapCC", "maskROI", "vectorAIF", "write.data.to.file")

        if(write.data.to.file==TRUE){
            save(dcemri.data, file=file.name)
            cat("...file saved as", file.name, "...")
            cat("\n")
            cat("...use the KAT() function to analyze data within this file.")
            cat("\n")
            cat("\n")
        }

        if(write.data.to.file==FALSE)
            return(dcemri.data)
    }
