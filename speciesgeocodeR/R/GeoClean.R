GeoClean <- function(x, isna = TRUE, isnumeric = TRUE, coordinatevalidity = TRUE, containszero = TRUE, zerozero = TRUE, zerozerothresh = 1,
                     latequallong = TRUE, GBIFhead = FALSE, countrycentroid = FALSE, contthresh = 0.5, capitalcoords = FALSE, capthresh = 0.5, 
                     countrycheck = FALSE, polygons, referencecountries = countryref, outp = c("summary", "detailed", "cleaned")) {
    
    dat <- x

    if ("lon" %in% names(dat)) {
        dat$XCOOR <- unlist(dat["lon"])
    }
    if ("lat" %in% names(dat)) {
        dat$YCOOR <- unlist(dat["lat"])
    }
    if ("longitude" %in% names(dat)) {
        dat$XCOOR <- unlist(dat["longitude"])
    }
    if ("latitude" %in% names(dat)) {
        dat$YCOOR <- unlist(dat["latitude"])
    }
    
    # old GBIF format
    if (dim(x)[2] == 25) {
        dat$identifier <- dat$species
        dat$country <- dat$ISO2
    }
    if (dim(x)[2] == 225) {
        dat$identifier <- dat$species
        dat$XCOOR <- dat$decimalLongitude
        dat$YCOOR <- dat$decimalLatitude
        dat$country <- dat$countryCode
        dat <- dat[c("identifier", "XCOOR", "YCOOR", "country")]
        dat <- data.frame(unlist(apply(dat, 2, function(x) gsub("^$|^ $", NA, x))))
    }
    
    if(dim(dat)[2] == 3){
      dat <- dat[c("identifier", "XCOOR", "YCOOR")]
      countrycheck <- FALSE
      warning("no country information found, countrycheck test skipped")
    }else{
    dat <- dat[c("identifier", "XCOOR", "YCOOR", "country")]
    }
    verb <- dat
    dat$clean <- T
    
    if (isna == T) {
        dat$clean[which(is.na(dat$XCOOR) | is.na(dat$YCOOR))] <- FALSE
        
        verb$isna <- T
        verb$isna[which(is.na(dat$XCOOR) | is.na(dat$YCOOR))] <- FALSE
    }
    
    if (isnumeric == T) {
        # is numeric
        dat$clean[which(suppressWarnings(is.na(as.numeric(as.character(dat$XCOOR)))))] <- FALSE
        dat$clean[which(suppressWarnings(is.na(as.numeric(as.character(dat$YCOOR)))))] <- FALSE
        
        verb$isnumeric <- TRUE
        verb$isnumeric[which(suppressWarnings(is.na(as.numeric(as.character(dat$XCOOR)))))] <- FALSE
        verb$isnumeric[which(suppressWarnings(is.na(as.numeric(as.character(dat$YCOOR)))))] <- FALSE
    }
    
    if (coordinatevalidity == T) {
        # -180 < long < 180
        dat$clean[which(suppressWarnings(as.numeric(as.character(dat$XCOOR))) > 180 | suppressWarnings(as.numeric(as.character(dat$XCOOR))) < 
            -180)] <- FALSE
        
        dat$clean[which(suppressWarnings(as.numeric(as.character(dat$YCOOR))) > 90 | suppressWarnings(as.numeric(as.character(dat$YCOOR))) < 
            -90)] <- FALSE
        
        verb$coordinatevalidity <- T
        verb$coordinatevalidity[which(suppressWarnings(as.numeric(as.character(dat$XCOOR))) > 180 | suppressWarnings(as.numeric(as.character(dat$XCOOR))) < 
            -180)] <- FALSE
        
        verb$coordinatevalidity[which(suppressWarnings(as.numeric(as.character(dat$YCOOR))) > 90 | suppressWarnings(as.numeric(as.character(dat$YCOOR))) < 
            -90)] <- FALSE
        
    }
    
    if (containszero == T) {
        # lat == 0 or long == 0
        dat$clean[which(suppressWarnings(as.numeric(as.character(dat$XCOOR))) == 0 | suppressWarnings(as.numeric(as.character(dat$YCOOR))) == 
            0)] <- FALSE
        verb$haszero <- TRUE
        verb$haszero[which(suppressWarnings(as.numeric(as.character(dat$XCOOR))) == 0 | suppressWarnings(as.numeric(as.character(dat$YCOOR))) == 
            0)] <- FALSE
        
    }
    if (zerozero == T) {

        loncap <- suppressWarnings(as.numeric(as.character(dat$XCOOR)) > (0 - zerozerothresh)) & 
                  suppressWarnings(as.numeric(as.character(dat$XCOOR))) <  (0 + zerozerothresh) 
        latcap <- suppressWarnings(as.numeric(as.character(dat$YCOOR)) > (0 - zerozerothresh)) & 
                  suppressWarnings(as.numeric(as.character(dat$YCOOR))) <  (0 + zerozerothresh) 

        dat$clean[which(loncap == T & latcap == T)] <- FALSE
        verb$zerozero <- TRUE
        verb$zerozero[which(loncap == T & latcap == T)] <- FALSE
    }
    
    if (latequallong == T) {
        # lat == long
        dat$clean[which(suppressWarnings(as.numeric(as.character(dat$XCOOR))) == suppressWarnings(as.numeric(as.character(dat$YCOOR))))] <- FALSE
        
        verb$latequallong <- TRUE
        verb$latequallong[which(suppressWarnings(as.numeric(as.character(dat$XCOOR))) == suppressWarnings(as.numeric(as.character(dat$YCOOR))))] <- FALSE
    }
    
    if (GBIFhead == T) {
        # degree around copenhagen
        loncop <- suppressWarnings(as.numeric(as.character(dat$XCOOR))) > 12.1 & suppressWarnings(as.numeric(as.character(dat$XCOOR))) < 
            12.8
        latcop <- suppressWarnings(as.numeric(as.character(dat$YCOOR))) > 55.5 & suppressWarnings(as.numeric(as.character(dat$YCOOR))) < 
            55.8
        dat$clean[which(loncop == T & latcop == T)] <- FALSE
        
        verb$GBIFhead <- T
        verb$GBIFhead[which(loncop == T & latcop == T)] <- FALSE
        
    }
    
    if (countrycentroid == T) {
        # 0.1 degree around country center
        countryref <- referencecountries
        conttest <- apply(dat, 1, function(x) .testcordcountr(x, countryref, contthresh))
        dat$clean[which(conttest == FALSE)] <- FALSE
        verb$countrycentroid <- conttest
    }
    if (capitalcoords == T) {
        # 0.1 degree around country capital #TESTTHIS
        countryref <- referencecountries
        captest <- apply(dat, 1, function(x) .testcordcap(x, countryref, capthresh))
        dat$clean[which(captest == FALSE)] <- FALSE
        verb$capitalcoordinates <- captest
    }
    if (countrycheck == T) {
        dat$XCOOR <- as.numeric(as.character(dat$XCOOR))
        dat$YCOOR <- as.numeric(as.character(dat$YCOOR))
        
        inp <- ReadPoints(dat[c("identifier", "XCOOR", "YCOOR")], polygons)
        
        if (all(nchar(as.character(dat$country)) <= 2, na.rm = T)) {
            contest <- SpGeoCodH(inp, areanames = "ISO2")
        }
        if (all(nchar(as.character(dat$country)) <= 3, na.rm = T) & !all(nchar(as.character(dat$country)) <= 2, na.rm = T)) {
            contest <- SpGeoCodH(inp, areanames = "ISO3")
        }
        if (!all(nchar(as.character(dat$country)) <= 3, na.rm = T)) {
            contest <- SpGeoCodH(inp, areanames = "NAME")
            warning("found country information with more than 3 letters; Country information should be ISO2 or ISO3")
        }
        
        if (!length(as.character(dat$country)) == length(as.character(contest$sample_table[, 2]))) {
            stop("coordinates include non-numerical or invalid elements; please check this before using the countrycheck argument")
        } else {
            verb$country.check <- as.character(dat$country) == as.character(contest$sample_table[, 2])
            dat$clean[which(verb$country.check == FALSE)] <- FALSE
            dat$clean[is.na(verb$country.check)] <- TRUE
            warning("unidentified country information test skipped and value set to TRUE without testing")
        }
    }
    if (outp[1] == "detailed") {
      verb$summary <- dat$clean
      rownames(verb) <- rownames(x)
        return(verb)
    }
    if(outp[1] == "summary"){
        return(dat$clean)
    }
    if(outp[1] == "cleaned"){
      verb$summary <- dat$clean
      rownames(verb) <- rownames(x)
      if(dim(x)[2] == 3){
        out <- verb[verb$summary == TRUE, 1:3]
      }else{ 
      out <- verb[verb$summary == TRUE, 1:4]
      }
      out$XCOOR <- as.numeric(as.character(out$XCOOR))
      out$YCOOR <- as.numeric(as.character(out$YCOOR))
      return(out)
    
    }

} 
