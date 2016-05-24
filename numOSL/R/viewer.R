##### 
viewer <- 
function(filename)  {
    UseMethod("viewer")
} #
### 2015.04.30.
viewer.default <- 
function(filename)  {
    ### Stop if not. 
    stopifnot(is.character(filename), 
              length(filename)==1L)
    ###
    binFile <- file(description=filename, open="rb")
    ###
    records <- list()
    iter <- 0L
    ### 
    ### Read each record iteratively.
    repeat  {
        ###  Data format version number (Version), Byte (2); ----------- USE.
        version <- readBin(binFile, what="integer", n=1L, size=2L, endian="little")
        if (length(version)==0L)  break 
        ### 
        iter <- iter + 1L
        ###
        ### Detailed description about BIN file data format  
        ### can be found from user mannual (version 3.22b)
        ### of the Analyst software, page 40, Appendix 1.
        ##################################################
        if (version!=6L)  {
              TotalLumType <- c("tl", "osl", "irsl", "m-ir", "m-vis", 
                                "tol", "pulsed", "rir", "rbr", "user")
              TotalDataType <- c("natural", "n+dose", "bleach", "bleach+dose", 
                                 "natural(bleach)", "n+dose(bleach)", "dose", "background")
              TotalLightSource <- c("none", "lamp", "ir-diodes",
                                    "calibraition-led", "blue-diodes")
              ### 
              ### Length of this record (Length), Small Integer (2);
              ### Length of previous record (Previous), Small Integer (2);
              ### Number of data points (NPoints), Small Integer (2). ----------- USE.
              pass <- readBin(binFile, what="integer", n=3L, size=2L, endian="little") 
              nPoints<- pass[3L]
              ###
              ### Luminescence type (LType), Byte (1).  ----------- USE.
              lumType <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
              ###
              ### Low (temperature, time, wavelength) (Low), Single (4); ----------- USE.
              ### High (temperature, time, wavelength) (High), Single (4); ----------- USE.
              ### Rate (heating rate, scan rate) (Rate), SIngle (4).
              pass <- readBin(binFile, what="double", n=3L, size=4L, endian="little")
              low <- pass[1L]
              high <- pass[2L]
              ###
              ###  Sample temperature (Temperature), Small Integer (2);
              ### X position of a single grain (XCoord), Small Integer (2);
              ### Y position of a single grain (YCOORD), Small Integer (2);
              ### TOL "delay" channels (Delay), Small Integer (2);
              ### TOL "on" channels (On), Small Integer (2);
              ### TOL "off" channels (Off), Small Integer (2).
              pass <- readBin(binFile, what="integer", n=6L, size=2L, endian="little")
              ###
              ### Carousel position (Position), Byte (1); ----------- USE.
              ### Run number (Run), Byte (1). ----------- USE.
              pass <- readBin(binFile, what="integer", n=2L, size=1L, endian="little")
              position <- pass[1L]
              runNumber <- pass[2L]
              ###
              ### Data collection time (Time), String (7); 
              ### Data collection date (Date), String (7). 
              pass <- readChar(binFile, nchars=14L, useBytes=TRUE)
              ###
              ### Sequence name (Sequence), String (9);
              ### User name (User), String (9).
              pass <- readChar(binFile, nchars=18L, useBytes=TRUE)
              ###
              ### Data type (Dtype), Byte (1). ----------- USE.
              dataType <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
              ###
              ### Irradiation time (IRR_Time), Single (4). ----------- USE.
              irrTime <- readBin(binFile, what="double", n=1L, size=4L, endian="little") 
              ###
              ### Irradiation type (alpha, beta or gamma) (IRR_Type), Byte (1); 
              ### Irradiation unit (Gy, Rads, secs, mins, hrs) (IRR_Unit), Byte (1). 
              pass <- readBin(binFile, what="integer", n=2L, size=1L, endian="little")
              ###
              ### Bleaching time (BL_Time), Single (4).    
              pass <- readBin(binFile, what="double", n=1L, size=4L, endian="little")
              ###
              ### Bleaching unit (mJ, J, secs, mins, hrs) (BL_Unit), Byte (1).    
              pass <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
              ###
              ### Annealing temperature (AN_Temp), Single (4);
              ### Annealing time (AN_Time), Single(4);
              ### Normalisation factor (1) (Norm1), Single (4);
              ### Normalisation factor (2) (Norm2), Single (4);
              ### Normalisation factor (3) (Norm3), Single (4);
              ### Background level (BG), Single (4).
              pass <- readBin(binFile, what="double", n=6L, size=4L, endian="little")
              ###
              ### Number of channels to shift data (Shift), Small Integer (2).    
              pass <- readBin(binFile, what="integer", n=1L, size=2L, endian="little")
              ###
              ### Sample name (Sample), String (21).
              pass <- readChar(binFile, nchars=21L, useBytes=TRUE) 
              ###
              ### Comment (Comment), String (81).
              pass <- readChar(binFile, nchars=81L, useBytes=TRUE) 
              ###
              ### Light Source (LightSource), Byte (1); ----------- USE.
              ### Ste Number (SET), Byte (1);
              ### Tag (TAG), Byte(1).
              pass <- readBin(binFile, what="integer", n=3L, size=1L, endian="little") 
              lightSource <- pass[1L]
              ###
              ### Grain number (Grain), Small Integer (2).    
              pass <- readBin(binFile, what="integer", n=1L, size=2L, endian="little")
              ###
              ### Optical stimulation power (LightPower), Single (4).    
              pass <- readBin(binFile, what="double", n=1L, size=4L, endian="little")  
              ###
              ### System ID (ystemID), Small Integer (2).
              pass <- readBin(binFile, what="integer", n=1L, size=2L, endian="little")
              ###
              ### Reserved (Reserved), Byte (54).    
              pass <- readBin(binFile, what="integer", n=54L, size=1L, endian="little")
              ###
              ### Data array of nPoints long integers (DPoints), Integer (4). ----------- USE.
              dataPoints <- readBin(binFile, what="integer", n=nPoints, size=4L, endian="little")
        } else {
              ### Detailed description about BINX file data format can be found 
              ### from User manual (February 2013) of Viewer software, page 38-39.
              #####################################################################
              TotalLumType <- c("tl", "osl", "irsl", "m-ir", "m-vis",  "tol", "trposl",
                                "rir", "rbr", "user", "posl", "sgosl", "rl", "xrf")
              TotalDataType <- c("natural", "n+dose", "bleach", "bleach+dose", 
                                 "natural(bleach)", "n+dose(bleach)", "dose", "background")
              TotalLightSource <- c("none", "lamp", "ir-diodes",
                                    "calibraition-led", "blue-diodes", "white-light", 
                                    "green-laser(sg)", "ir-laser(sg)")
              ###
              ### Length of this record (Length), Long Integer (4);
              ### Length of previous record (Previous), Long Integer (4);
              ### Number of data points (NPoints), Long Integer (4). ----------- USE.
              pass <- readBin(binFile, what="integer", n=3L, size=4L, endian="little") 
              nPoints <- pass[3L]
              ###
              ### Run number (Run), Small Integer (2); ----------- USE.
              ### Set number (Set), Small Integer (2);
              ### Carousel position (Position), Small Integer (2); ----------- USE.
              ### Grain number (GrainNumber), Small Integer (2);
              ### Curve number (for multiple curve operations) (CurveNo), Small Integer (2);
              ###  X position of a single grain (Xcoord), Small Integer (2);
              ###  Y position of a single grain (Ycoord), Small Integer (2). 
              pass <- readBin(binFile, what="integer", n=7L, size=2L, endian="little")     
              runNumber <- pass[1L]
              position <- pass[3L] 
              ###
              ### Sample name (Sample), String (21).
              pass <- readChar(binFile, nchars=21L, useBytes=TRUE) 
              ###
              ### Comment (comment), String (81).
              pass <- readChar(binFile, nchars=81L, useBytes=TRUE) 
              ###
              ### System ID (SystemID), Small Integer (2). 
              pass <- readBin(binFile, what="integer", n=1, size=2L, endian="little") 
              ###
              ### File name (.SEC, BINX ect) (FName), String (101).
              pass <- readChar(binFile, nchars=101L, useBytes=TRUE)
              ###
              ### User name (User), String (31).
              pass <- readChar(binFile, nchars=31L, useBytes=TRUE)
              ###
              ### Data collection time (hh-mm-ss) (Time), String (7);
              ### Data collection date (dd-mm-yy) (Date), String (7).
              pass <- readChar(binFile, nchars=14L, useBytes=TRUE)
              ###
              ### Data type (DType), Byte (1). ----------- USE.
              dataType <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
              ###
              ### Bleaching time (BL_Time), Single (4).     
              pass <- readBin(binFile, what="double", n=1L, size=4L, endian="little")
              ###
              ### Bleaching unit (mJ, J, secs, mins, hrs) (BL_Unit), Byte (1).    
              pass <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
              ###
              ### Normalisation factor (1) (Norm1), Single (4);
              ### Normalisation factor (2) (Norm2), Single (4);
              ### Normalisation factor (3) (Norm3), Single (4);
              ### Background level (BG), single (4).
              pass <- readBin(binFile, what="double", n=4, size=4L, endian="little") 
              ### 
              ### Number of channels to shift data (Shift), Small Integer (2).    
              pass <- readBin(binFile, what="integer", n=1L, size=2L, endian="little")
              ###
              ### Tag (Tag), Byte (1).
              pass <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
              ###   
              ### Reserved for internal use (??Reserved),  ?? (20).
              pass <- readBin(binFile, what="raw", n=20L, size=1L, endian="little")
              ###
              ### Luminescence type (LType), Byte (1); ----------- USE.
              ### Light Source (LightSource), Byte (1).  ----------- USE.
              pass <- readBin(binFile, what="integer", n=2L, size=1L, endian="little")
              lumType <- pass[1L]
              lightSource <- pass[2L]
              ###
              ### Optical Stimulation Power (LightPower), Single (4); 
              ### Low (temperature, time, wavelength) (Low), Single (4); ----------- USE.
              ### High (temperature, time, wavelength) (High), Single (4); ----------- USE.
              ### Rate (heating rate, scan rate) (Rate), Single (4).
              pass <- readBin(binFile, what="double", n=4L, size=4L, endian="little")
              low <- pass[2L]
              high <- pass[3L]
              ###
              ### Sample temperature (Temperature), Small Integer (2);
              ### Measurement temperature (MeasTemp),Small Integer (2).
              pass <- readBin(binFile, what="integer", n=2L, size=2L, endian="little")
              ###
              ### Preheating temperature (An_Temp), Single (4);
              ### Preheating time (An_Time), Single (4).
              pass <- readBin(binFile, what="double", n=2L, size=4L, endian="little")
              ###
              ### TOL "delay" channels (Delay), Small Integer (2);
              ### TOL "on" channels (On), Small Integer (2);
              ### TOL "off" channels (off), Small Integer (2). 
              pass <- readBin(binFile, what="integer", n=3L, size=2L, endian="little")
              ###
              ### Irradiation time (IRR_Time), Single (4). ----------- USE.
              irrTime <- readBin(binFile, what="double", n=1L, size=4L, endian="little")
              ###
              ### Irradiation type (alpha, beta or gamma) (IRR_Type), Byte (1).
              pass <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
              ### 
              ### Irradiation dose rate (Gy/s) (IRR_DoseRate), Single (4);
              ### Irradiation dose rate error (Gy/s) (DoseRateErr), Single (4).
              pass <- readBin(binFile, what="double", n=2L, size=4L, endian="little")
              ###
              ### Time since last irradiation (s) (TimeSinceIrr), Long Integer (4);
              pass <- readBin(binFile, what="integer", n=1L, size=4L, endian="little")
              ###     
              ### Time unit (time tick) for pulse parameters (s) (TimeTick), Single (4).
              pass <- readBin(binFile, what="double", n=1L, size=4L, endian="little")
              ###
              ### On-time for pulsed stimulation (in time ticks) (OnTime), Long Integer (4);
              ### Stimulation period (on+off time in time ticks) (StimPeriod), Long Integer (4).
              pass <- readBin(binFile, what="integer", n=2L, size=4L, endian="little")
              ###  
              ### PMT signal gating enabled (GateEnabled), Byte (1).
              pass <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
              ###
              ### Start of gating (in time ticks from start of on pulse) (GateStart), Long Integer (4);
              ### End of gating (in time ticks from start of on pulse) (GateEnd), Long Integer (4).
              pass <- readBin(binFile, what="integer", n=2L, size=4L, endian="little")
              ###
              ### Photon Timer enabled (PTenabled), Byte (1);
              ### PMT dead time correction enabled (DTenabled), Byte (1).
              pass <- readBin(binFile, what="integer", n=2L, size=1L, endian="little")
              ### 
              ### PMT dead time (s) (DeadTime), Single (4);
              ### Stimulation power corresponding to 100% (mw/cm2) (MaxLPower), Single (4);
              ### XRF acquisition time (s) (XrfAcqTime), Single (4);
              ### XRF X-ray high voltage (V) (XrfHV), Single (4). 
              pass <- readBin(binFile, what="double", n=4L, size=4L, endian="little")
              ###
              ### XRF X-ray current (uA) (XrfCurr), Long Integer (4).
              pass <- readBin(binFile, what="integer", n=1L, size=4L, endian="little")
              ###
              ### XRF dead time fraction (XrfDeadTimeF), Single (4).
              pass <- readBin(binFile, what="double", n=1L, size=4L, endian="little")
              ###
              ### Reserved for internal use (??Reserved1),   Byte (24). 
              pass <- readBin(binFile, what="raw", n=24L, size=1L, endian="little")
              ###
              ### Data array of nPoints Long Integers (DPoints), Long Integer (4).----------- USE.
              dataPoints <- readBin(binFile, what="integer", n=nPoints, size=4L, endian="little")
       } # end if.
       ###
       attr(dataPoints, "position") <- position 
       attr(dataPoints, "runNumber") <- runNumber
       attr(dataPoints, "nPoints") <- nPoints
       attr(dataPoints, "low") <- low
       attr(dataPoints, "high") <- high
       attr(dataPoints, "irrTime") <- round(irrTime, digits=3L)
       attr(dataPoints, "lumType") <- TotalLumType[lumType+1L]
       attr(dataPoints, "dataType") <- TotalDataType[dataType+1L]     
       attr(dataPoints, "lightSource") <- TotalLightSource[lightSource+1L]
       ###
       records[[iter]] <- dataPoints
    } # end repeat.
    ### 
    close(binFile)
    ###
    tab <- data.frame("position"=sapply(X=records, FUN=attr, "position"),
                      "runNumber"=sapply(X=records, FUN=attr, "runNumber"),
                      "nPoints"=sapply(X=records, FUN=attr, "nPoints"),
                      "low"=sapply(X=records, FUN=attr, "low"),
                      "high"=sapply(X=records, FUN=attr, "high"),
                      "irrTime"=sapply(X=records, FUN=attr, "irrTime"),
                      "lumType"=sapply(X=records, FUN=attr, "lumType"),
                      "dataType"=sapply(X=records, FUN=attr, "dataType"),
                      "lightSource"=sapply(X=records, FUN=attr, "lightSource"))
    ###
    cat(paste("Number of records: ", length(records), "\n\n", sep=""))
    ###
    cat("Luminescence type(s):\n")
    print(levels(factor(sapply(records, attr, "lumType"))))
    cat("\n")
    ###
    cat("Carousel position(s):\n")
    print(as.numeric(levels(factor(sapply(records, attr, "position")))))
    cat("\n")
    ###
    cat("Data Type(s):\n")
    print(levels(factor(sapply(records, attr, "dataType"))))
    cat("\n")
    ###
    cat("Levels of irradiation time:\n")
    print(as.numeric(levels(factor(sapply(records, attr, "irrTime")))))
    cat("\n")
    ###
    cat("Light source(s):\n")
    print(levels(factor(sapply(records, attr, "lightSource"))))
    cat("\n")
    ### 
    output <- list("records"=records, "tab"=tab)
    class(output) <- "binfile"
    invisible(output)
} # end function viewer.
#####
