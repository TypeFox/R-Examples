## conversion of absolute size units
getConvFac <-
function(conversion="m2cm") {
    ## check how the conversion factor is indicated
    convFac <- if(is.character(conversion)) {
        idxMM2MM <- conversion ==     "mm2mm"
        idxMM2CM <- conversion ==     "mm2cm"
        idxMM2M  <- conversion ==     "mm2m"
        idxMM2KM <- conversion ==     "mm2km"
        idxMM2IN <- conversion %in% c("mm2in", "mm2inch")
        idxMM2FT <- conversion %in% c("mm2ft", "mm2feet")
        idxMM2YD <- conversion %in% c("mm2yd", "mm2yard")

        idxCM2MM <- conversion ==     "cm2mm"
        idxCM2CM <- conversion ==     "cm2cm"
        idxCM2M  <- conversion ==     "cm2m"
        idxCM2KM <- conversion ==     "cm2km"
        idxCM2IN <- conversion %in% c("cm2in", "cm2inch")
        idxCM2FT <- conversion %in% c("cm2ft", "cm2feet")
        idxCM2YD <- conversion %in% c("cm2yd", "cm2yard")

        idxM2MM  <- conversion ==     "m2mm"
        idxM2CM  <- conversion ==     "m2cm"
        idxM2M   <- conversion ==     "m2m"
        idxM2KM  <- conversion ==     "m2km"
        idxM2IN  <- conversion %in% c("m2in", "m2inch")
        idxM2FT  <- conversion %in% c("m2ft", "m2feet")
        idxM2YD  <- conversion %in% c("m2yd", "m2yard")

        idxKM2MM <- conversion ==     "km2mm"
        idxKM2CM <- conversion ==     "km2cm"
        idxKM2M  <- conversion ==     "km2m"
        idxKM2KM <- conversion ==     "km2km"
        idxKM2IN <- conversion %in% c("km2in", "km2inch")
        idxKM2FT <- conversion %in% c("km2ft", "km2feet")
        idxKM2YD <- conversion %in% c("km2yd", "km2yard")

        idxIN2MM <- conversion %in% c("in2mm", "inch2mm")
        idxIN2CM <- conversion %in% c("in2cm", "inch2cm")
        idxIN2M  <- conversion %in% c("in2m",  "inch2m")
        idxIN2KM <- conversion %in% c("in2km", "inch2km")
        idxIN2IN <- conversion %in% c("in2in", "in2inch", "inch2in", "inch2inch")
        idxIN2FT <- conversion %in% c("in2ft", "in2feet", "inch2ft", "inch2feet")
        idxIN2YD <- conversion %in% c("in2yd", "inch2yard")

        idxYD2MM <- conversion %in% c("yd2mm", "yard2mm")
        idxYD2CM <- conversion %in% c("yd2cm", "yard2cm")
        idxYD2M  <- conversion %in% c("yd2m",  "yard2m")
        idxYD2KM <- conversion %in% c("yd2km", "yard2km")
        idxYD2IN <- conversion %in% c("yd2in", "yard2inch")
        idxYD2FT <- conversion %in% c("yd2ft", "yd2feet", "yard2ft", "yard2feet")
        idxYD2YD <- conversion %in% c("yd2yd", "yd2yard", "yard2yd", "yard2yard")

        idxFT2MM <- conversion %in% c("ft2mm", "feet2mm")
        idxFT2CM <- conversion %in% c("ft2cm", "feet2cm")
        idxFT2M  <- conversion %in% c("ft2m",  "feet2m")
        idxFT2KM <- conversion %in% c("ft2km", "feet2km")
        idxFT2IN <- conversion %in% c("ft2in", "ft2inch", "feet2in", "feet2inch")
        idxFT2FT <- conversion %in% c("ft2ft", "feet2feet")
        idxFT2YD <- conversion %in% c("ft2yd", "feet2yd", "ft2yard", "feet2yard")

        ## did we catch all requested conversion types?
        idxAll <- idxMM2MM | idxMM2CM | idxMM2M | idxMM2KM | idxMM2IN | idxMM2FT | idxMM2YD |
                  idxCM2MM | idxCM2CM | idxCM2M | idxCM2KM | idxCM2IN | idxCM2FT | idxCM2YD |
                  idxM2MM  | idxM2CM  | idxM2M  | idxM2KM  | idxM2IN  | idxM2FT  | idxM2YD  |
                  idxKM2MM | idxKM2CM | idxKM2M | idxKM2KM | idxKM2IN | idxKM2FT | idxKM2YD |
                  idxIN2MM | idxIN2CM | idxIN2M | idxIN2KM | idxIN2IN | idxIN2FT | idxIN2YD |
                  idxFT2MM | idxFT2CM | idxFT2M | idxFT2KM | idxFT2IN | idxFT2FT | idxFT2YD |
                  idxYD2MM | idxYD2CM | idxYD2M | idxYD2KM | idxYD2IN | idxYD2FT | idxYD2YD
        if(!all(idxAll)) {
            warning(c('Conversion type(s) "', paste(conversion[!idxAll], collapse=", "),
                      '" not found - conversion factor set to 1'))
        }

        res <- rep(1, length(conversion))

        ## conversion factors for length units
        res[idxMM2MM] <-       1
        res[idxMM2CM] <-       1/10
        res[idxMM2M]  <-       1/1000
        res[idxMM2KM] <-       1/1000000
        res[idxMM2IN] <-       1/25.4
        res[idxMM2FT] <-       1/304.8
        res[idxMM2YD] <-       1/914.4

        res[idxCM2MM] <-      10
        res[idxCM2CM] <-       1
        res[idxCM2M]  <-       1/100
        res[idxCM2KM] <-       1/100000
        res[idxCM2IN] <-       1/2.54
        res[idxCM2FT] <-       1/30.48
        res[idxCM2YD] <-       1/91.44

        res[ idxM2MM] <-    1000
        res[ idxM2CM] <-     100
        res[ idxM2M]  <-       1
        res[ idxM2KM] <-       1/1000
        res[ idxM2IN] <-       1/0.0254
        res[ idxM2FT] <-       1/0.3048
        res[ idxM2YD] <-       1/0.9144

        res[idxKM2MM] <- 1000000
        res[idxKM2CM] <-  100000
        res[idxKM2M]  <-    1000
        res[idxKM2KM] <-       1/1
        res[idxKM2IN] <-       1/0.0000254
        res[idxKM2FT] <-       1/0.0003048
        res[idxKM2YD] <-       1/0.0009144

        res[idxIN2MM] <-      25.4
        res[idxIN2CM] <-       2.54
        res[idxIN2M]  <-       0.0254
        res[idxIN2KM] <-       0.0000254
        res[idxIN2IN] <-       1
        res[idxIN2FT] <-       1/12
        res[idxIN2YD] <-       1/36

        res[idxFT2MM] <-     304.8
        res[idxFT2CM] <-      30.48
        res[idxFT2M]  <-       0.3048
        res[idxFT2KM] <-       0.0003048
        res[idxFT2IN] <-      12
        res[idxFT2FT] <-       1
        res[idxFT2YD] <-       1/3

        res[idxYD2MM] <-     914.4
        res[idxYD2CM] <-      91.44
        res[idxYD2M]  <-       0.9144
        res[idxYD2KM] <-       0.0009144
        res[idxYD2IN] <-      36
        res[idxYD2FT] <-       3
        res[idxYD2YD] <-       1

        res
    } else if(is.numeric(conversion)) {  # factor is given directly
        if(any(conversion <= 0)) {
            stop("Conversion factor must be > 0")
        }
        conversion
    } else {
        stop("Conversion must be a character constant or numeric factor")
    }                                # if(is.character(conversion))

    return(convFac)
}

## determine unit of (x,y)-coords from conversion string
getUnits <-
function(x="m2cm", first=TRUE) {
    if(!is.character(x)) {
        warning("Unit not recognized - input must have form like m2cm")
        return(" ")
    }

    units    <- strsplit(x, "2")                  # first and second part of string
    unitLens <- vapply(units, length, integer(1)) # count parts
    if(!all(unitLens == 2L)) {                    # check that there are two parts
        warning("Unit not recognized - input must have form like m2cm")
        return("")
    }

    knownUnits <- c("m", "cm", "mm", "yd", "yard", "ft", "foot", "feet", "in", "inch")
    isKnown    <- vapply(units, function(x) { all(x %in% knownUnits) }, logical(1))
    if(!all(isKnown)) {
        warning(c("Unit not recognized - needs to be one of\n",
                  paste(knownUnits, collapse=" ")))
        return("")
    }

    ## replace feet with ft, yard with yd, inch with in
    replaceUnit <- function(x) {
        x <- gsub("yard", "yd", x)
        x <- gsub("foot", "ft", x)
        x <- gsub("feet", "ft", x)
        gsub("inch", "in", x)
    }

    units <- lapply(units, replaceUnit)
    if(first) {
        vapply(units, head, FUN.VALUE=character(1), n=1)
    } else {
        vapply(units, tail, FUN.VALUE=character(1), n=1)
    }
}
