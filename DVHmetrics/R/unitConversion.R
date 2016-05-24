## conversion of dose and volume units
getConvFac <-
function(conversion="CGY2GY") {
    conversion <- toupper(removeWS(conversion))

    ## check how the conversion factor is indicated
    idxGY2GY   <- conversion %in% c("GY2GY")
    idxCGY2GY  <- conversion %in% c("CGY2GY")
    idxGY2CGY  <- conversion %in% c("GY2CGY")
    idxCGY2CGY <- conversion %in% c("CGY2CGY")
    idxCC2CC   <- conversion %in% c("CC2CC")

    ## did we catch all requested conversion types?
    idxAll <- idxGY2GY | idxCGY2GY | idxGY2CGY | idxCGY2CGY | idxCC2CC
    if(!all(idxAll)) {
        warning(c('Conversion type(s) "', paste(conversion[!idxAll], collapse=", "),
                  '" not found - conversion factor set to NA'))
    }

    convFac <- rep(NA_real_, length(conversion))

    ## conversion factors for dose units
    convFac[idxGY2GY]   <- 1
    convFac[idxCGY2GY]  <- 1/100
    convFac[idxGY2CGY]  <- 100
    convFac[idxCGY2CGY] <- 1
    convFac[idxCC2CC]   <- 1

    return(convFac)
}

## determine unit from conversion string
getUnits <-
function(x="CGY2GY", first=TRUE) {
    if(!is.character(x)) {
        warning("Unit not recognized - input must have form like CGY2GY")
        return(" ")
    }

    units    <- strsplit(x, "2")         # first and second part of string
    unitLens <- sapply(units, length)    # count parts
    if(!all(unitLens == 2)) {            # check that there are two parts
        warning("Unit not recognized - input must have form like CGY2GY")
        return(NA_character_)
    }

    knownUnits <- c("CGY", "GY", "CC")
    isKnown    <- sapply(units, function(x) { all(x %in% knownUnits) })
    if(!all(isKnown)) {
        warning(c("Unit not recognized - needs to be one of\n",
                  paste(knownUnits, collapse=" ")))
        return("")
    }

    if(first) {
        sapply(units, head, n=1)
    } else {
        sapply(units, tail, n=1)
    }
}
