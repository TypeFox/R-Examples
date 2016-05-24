GribInfo <- function(grib.file, file.type = "grib2") {
    #This function returns information about what's contained in a grib file.
    #INPUTS
    #    GRIB.FILE - Path and name of file to examine
    #    FILE.TYPE = Whether it's a grib2 file ("grib2") or a grib file ("grib1")
    #OUTPUTS
    #    GRIB.INFO - What the grib file contains
    #        $INVENTORY - Information on variables, levels, and forecasts
    #        $GRID - Information on the model grid, only available in grib2 
 
    if(file.type == "grib2") {
        op <- options("warn")
        options(warn = -1)
        test <- tryCatch(system('wgrib2', intern = TRUE))
        options(op)
        if(attr(test, "status") != 8) {
            stop("wgrib2 does not appear to be installed, or it is not on the PATH variable.
                You can find wgrib2 here: http://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/.
                If the binaries don't work, try compiling from source.")
        }
        inv <- system(paste0("wgrib2 ", grib.file, " -inv -"), intern = TRUE)
        grid <- system(paste0("wgrib2 ", grib.file, " -grid"), intern = TRUE) 
    } else if (file.type == "grib1") {
         op <- options("warn")
         options(warn = -1)
         test <- tryCatch(system('wgrib', intern = TRUE))
         options(op)
          if(attr(test, "status") != 8) {
              stop("wgrib does not appear to be installed, or it is not on the PATH variable.
                  You can find wgrib here: http://www.cpc.ncep.noaa.gov/products/wesley/wgrib.html.")
           }
        inv <- system(paste0("wgrib ", grib.file), " -s", intern = TRUE)
        grid <- NULL
    } else {
        stop(paste0("Did not recognise file type ", file.type, ".  Please use \"grib2\" or \"grib.\""))
    }
    return(list(inventory = inv, grid = grid))
}

ReadGrib <- function(file.name, levels, variables, domain = NULL, file.type = "grib2", missing.data = NULL) {
    #This is a function to read forecast data from a Grib file
    #INPUTS
    #    FILE.NAME - Grib file name
    #    VARIABLES - data to extract
    #    LEVELS - which levels to extract data from
    #    DOMAIN - Region to extract data from, in c(LEFT LON, RIGHT LON, TOP LAT, BOTTOM LAT), west negative
    #    FILE.TYPE - whether this is a grib1 or a grib2 file
    #        If grib1, you must have the wgrib program installed
    #        If grib2, you must have the wgrib2 program installed
    #    MISSING.DATA - Replace missing data in grib archive with this value.
    #        If NULL, leave the data out.  Only works with wgrib2. Default NULL.
    #        See Trick 19 here: http://www.ftp.cpc.ncep.noaa.gov/wd51we/wgrib2/tricks.wgrib2
    #OUTPUTS
    #    MODEL.DATA - the grib model as an array, with columns for the model run date (when the model was run)
    #       the forecast (when the model was for), the variable (what kind of data), the level (where in the atmosphere or the Earth, vertically)
    #       the longitude, the latitude, and the value of the variable.

    #Get specified data from grib file

    if(file.type == "grib2") {
        op <- options("warn")
        options(warn = -1)
        test <- tryCatch(system('wgrib2', intern = TRUE))
        options(op)
        if(attr(test, "status") != 8) {
            stop("wgrib2 does not appear to be installed, or it is not on the PATH variable.
                You can find wgrib2 here: http://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/.
                If the binaries don't work, try compiling from source.")
        }
        match.str <- ' -match "('
        for(var in variables) {
            match.str <- paste(match.str, var, "|", sep = "")
        }
    
        match.str.lst <- strsplit(match.str, split = "")[[1]]
        match.str <- paste(match.str.lst[1:(length(match.str.lst) - 1)], collapse = "")
    
        if(length(levels) > 0 & !is.null(levels)) {
            match.str <- paste(match.str, "):(", sep = "")
            for(lvl in levels) {
                match.str <- paste(match.str, lvl, "|", sep = "")
            }
        } else {
            match.str <- paste0(match.str, ")")
       }
    
        match.str.lst <- strsplit(match.str, split = "")[[1]]
        match.str <- paste(match.str, '"', sep = "")
        match.str <- paste(match.str.lst[1:(length(match.str.lst) - 1)], collapse = "")
        match.str <- paste(match.str, ")\"", sep = "")
   
        if(!is.null(missing.data) & !is.numeric(missing.data)) {
            warning(paste("Your value", missing.data, " for missing data does not appear to be a number!"))
        }
        if(!(is.null(missing.data))) {
            missing.data.str <- paste0(" -rpn \"sto_1:", missing.data, ":rcl_1:merge\"")
        } else {
            missing.data.str <- ""
        }

        #Write out a grib file with a smaller domain, then read it in
        if(!is.null(domain)) {
           if(!length(domain) == 4 | any(!is.numeric(domain))) {
              stop("Input \"domain\" is the wrong length and/or consists of something other than numbers.
                  It should be a 4 element vector: c(LEFT LON, RIGHT LON, TOP LAT, BOTTOM LAT)")
           } else {
               wg2.pre <- paste0('wgrib2 ',
                   file.name,
                   " -small_grib ", 
                   domain[1], ":", domain[2], " ", domain[4], ":", domain[3],
                   " - | wgrib2 - ")
           }
        } else {
           wg2.pre <- paste0('wgrib2 ',  file.name)
        }
        
        wg2.str <- paste(wg2.pre,
            ' -inv my.inv',
            missing.data.str,
            ' -csv - -no_header', 
            match.str, sep = "")
        
        #Get the data from the grib file in CSV format
        if(Sys.info()[["sysname"]] == "Windows") {
            csv.str <- shell(wg2.str, intern = TRUE)        
        } else {
            csv.str <- system(wg2.str, intern = TRUE)
        } 

        #HERE IS THE EXTRACTION
        model.data.vector <- strsplit(paste(gsub("\"", "", csv.str), collapse = ","), split = ",")[[1]]
        chunk.inds <- seq(1, length(model.data.vector) - 6, by = 7)
        model.data <- list(model.run.date = model.data.vector[chunk.inds],
            forecast.date = model.data.vector[chunk.inds + 1],
            variables = model.data.vector[chunk.inds + 2],
            levels = model.data.vector[chunk.inds + 3],
            lon = as.numeric(model.data.vector[chunk.inds + 4]),
            lat = as.numeric(model.data.vector[chunk.inds + 5]),
            value = model.data.vector[chunk.inds + 6],
            meta.data = "None - this field is used for grib1 files",
            grib.type = file.type
            )
      } else if (file.type == "grib1") {
         op <- options("warn")
         options(warn = -1)
         test <- tryCatch(system('wgrib', intern = TRUE))
         options(op)
          if(attr(test, "status") != 8) {
              stop("wgrib does not appear to be installed, or it is not on the PATH variable.
                  You can find wgrib here: http://www.cpc.ncep.noaa.gov/products/wesley/wgrib.html.
                  It is also available as an Ubuntu package.")
           }
           
           # wgrib -s fcst.grb1 | grep ":TMP:1000 mb:" | wgrib -i -text fcst.grb1 -o asciifile.txt
           #This is inelegant - but I think it will work
            
           model.data <- list(meta.data = NULL, value = NULL, variables = NULL, levels = NULL)
           c <- 1
           for(var in variables) {
               for(lvl in levels) {
                   wg.str <- paste0("wgrib -s ", file.name, " | grep \":", 
                       var, ":", lvl, ":\" | wgrib -V -i -text ", file.name, " -o tmp.txt")
                   #The meta.data variable contains info on the lat/lon grid
                   model.data$meta.data[[c]] <- system(wg.str, ignore.stderr = TRUE)
                   model.data$value[[c]] <- scan("tmp.txt", skip = 1, quiet = TRUE) 
                   model.data$variables[c] <- var
                   model.data$levels[c] <- lvl 
                   c <- c + 1
               }
           }
           model.data$grib.type <- file.type    
    }
    return(model.data)
}
