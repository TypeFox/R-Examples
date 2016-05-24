MODISSummaries <-
function(LoadDat, FileSep = NULL, Dir = ".", Product, Bands, ValidRange, NoDataFill, ScaleFactor, StartDate = FALSE,
         QualityScreen = FALSE, QualityBand = NULL, QualityThreshold = NULL, Mean = TRUE, SD = TRUE, Min = TRUE, Max = TRUE,
         Yield = FALSE, Interpolate = FALSE, InterpolateN = NULL, DiagnosticPlot = FALSE) 
{ 
    # DEFINE
    NUM_METADATA_COLS <- 10
        
    if(Dir == '.') cat('MODIS data files from ', getwd(),
                       ' will be summarised.\nSummary files will be written to the same directory.\n', sep = '') 
    if(Dir != '.') cat('MODIS data files from ', Dir,
                       ' will be summarised.\nSummary files will be written to the same directory.\n', sep = '')
    
    # Load input time-series data file; external data file, or an R object.
    if(!is.object(LoadDat) & !is.character(LoadDat)) stop("LoadDat must be an object in R or a file path character string.")
    if(is.object(LoadDat)) details <- data.frame(LoadDat)
    if(is.character(LoadDat)){
      if(!file.exists(LoadDat)) stop("Character string input for LoadDat does not resemble an existing file path.")
      if(is.null(FileSep)) stop("To load a file as input, you must also specify its delimiter (FileSep).")
      details <- read.delim(LoadDat, sep = FileSep) 
    }
    
    #####
    if(!file.exists(Dir)) stop("Character string input for Dir argument does not resemble an existing file path.")
    
    # Check valid inputs for the quality screening process.
    if(!is.logical(QualityScreen)) stop("QualityScreen argument should be a logical input. See help ?MODISSummaries.")
    if(QualityScreen){
      if(is.null(QualityBand) | is.null(QualityThreshold)) stop("QualityBand and QualityThreshold not specified.")
    }
    
    product.bands <- try(GetBands(Product), silent = TRUE)
    if(class(product.bands) != "try-error"){
      # Check Band and QualityBand belong to Product.
      if(!all(Bands %in% product.bands)) stop(paste("Band input does not match with", Product, "product.", sep = " "))
      if(QualityScreen){
        if(Product == "MCD43A4"){
          if(QualityBand != "BRDF_Albedo_Band_Quality") stop("QualityBand input is not QA data for MCD43A4 product.")
        } else {
          if(!any(product.bands == QualityBand)) stop(paste("QualityBand is not QA data for", Product, "product.", sep = " "))
        }
      }
    } else {
      cat("MODIS server temporarily overloaded so user input checking skipped.")
    }
    
    # NoDataFill should be one integer.
    if(length(NoDataFill) != 1) stop("NoDataFill input must be one integer.")
    if(!is.numeric(NoDataFill)) stop("NoDataFill should be numeric class. One integer.")
    if(abs(NoDataFill - round(NoDataFill)) > .Machine$double.eps^0.5) stop("NoDataFill input must be one integer.")
    
    # ValidRange should be a numeric vector, length 2.
    if(length(ValidRange) != 2) stop("ValidRange input must be a numeric vector - an upper and lower bound.")
    if(!is.numeric(ValidRange)) stop("ValidRange should be numeric class.")
    
    # ScaleFactor should be numeric, length 1.
    if(length(ScaleFactor) != 1) stop("ScaleFactor input must be one numeric element.")
    if(!is.numeric(ScaleFactor)) stop("ScaleFactor should be numeric class.")
    
    # Year or posixt date format?
    if(all(is.na(details$end.date))) stop("Time-series end.date variable in LoadDat is empty.")
    
    Year <- FALSE
    POSIXt <- FALSE
    
    posix.compatible <- try(as.POSIXlt(details$end.date), silent=TRUE)
    
    if(any(class(details$end.date) == "POSIXt") | all(class(posix.compatible) != "try-error")) POSIXt <- TRUE
    if(all(is.numeric(details$end.date) & nchar(as.character(details$end.date)) == 4) & 
         any(class(posix.compatible) == "try-error")) Year <- TRUE
    
    if(!Year & !POSIXt) stop("Date informDate information in LoadDat is not recognised as years or as POSIXt format.")
    if(Year & POSIXt) stop("Date information in LoadDat is recognised as both year and POSIXt formats.")
    
    if(StartDate){
      if(!any(names(details) == "start.date") | all(is.na(details$start.date))){
        stop("StartDate == TRUE, but no start.date field found in LoadDat.")
      }
    }
    #####
    
    # Date and time stamp
    date <- as.POSIXlt(Sys.time())
    file.date <- paste(as.Date(date),
                       paste(paste0("h",date$hour), paste0("m",date$min), paste0("s",round(date$sec,digits=0)), sep = "-"),
                       sep = "_")
    
    # Get a list of all downloaded subset (.asc) files in the data directory.
    file.list <- list.files(path = Dir, pattern = paste(Product, ".*asc$", sep = ""))
    if(length(file.list) == 0) stop("Found no MODIS data files in Dir that match the request.")
    
    size.test <- sapply(file.path(Dir, file.list), function(x) ncol(read.csv(x)[1, ]) - NUM_METADATA_COLS)
    if(!all(size.test == size.test[1])) stop("The number of pixels (Size) in subsets identified are not all the same.")
    # size.test checked all subsets are compatible for processing to one summary file. Can now just use size.test[1]
    num.pixels <- unname(size.test[1])
    
    # Extract IDs for ASCII files so they can be included in summary output; ncol = length(file.list), nrow = size.test.
    where.id <- regexpr("___", file.list)

    id <- matrix(unname(mapply(function(x, y, z) rep(substr(x, 1, y - 1), z), x = file.list, y = where.id, z = size.test,
                               SIMPLIFY = "array")), nrow = size.test, ncol = length(file.list))

    # Time-series analysis for each time-series (ASCII file) consecutively.
    band.data.site <- lapply((size.test *  length(Bands)), function(x) matrix(nrow = x, ncol = 13))
    
    for(counter in 1:length(file.list)){
      
      cat("Processing file ", counter, " of ", length(file.list), "...\n", sep = "")
      
      ##### Load selected .asc file into a data frame, name columns and tell user what's being processed.
      ds <- read.csv(file.path(Dir, file.list[counter]), header = FALSE, as.is = TRUE)
      names(ds) <- c("nrow", "ncol", "xll", "yll", "pixelsize", "row.id", "product.code", "MODIS.acq.date",
                     "where", "MODIS.proc.date", 1:(ncol(ds) - NUM_METADATA_COLS))
      
      # Extract year and day from the metadata and make POSIXlt dates (YYYY-MM-DD), ready for time-series analysis.
      year <- as.numeric(substr(ds$MODIS.acq.date, 2, 5))
      day <- as.numeric(substr(ds$MODIS.acq.date, 6, 8))
      date <- strptime(paste(year, "-", day, sep = ""), "%Y-%j")
      ds <- cbind(ds[,1:NUM_METADATA_COLS], date, ds[,(NUM_METADATA_COLS+1):ncol(ds)]) 
      w.ds.dat <- which(names(ds) == "date") + 1
      #####
      
      ##### Section that uses the files metadata strings [ ,1:5] for each time-series to extract necessary information.
      wherelong <- regexpr("Lon", ds$where[1])
      whereSamp <- regexpr("Samp", ds$where[1])
      lat <- as.numeric(substr(ds$where[1], 4, wherelong - 1))
      long <- as.numeric(substr(ds$where[1], wherelong + 3, whereSamp - 1))
      
      # Check that all bands listed are in the ASCII files.
      Band.check <- sapply(Bands, function(x) any(grepl(x, ds$row.id)))
      if(!all(Band.check)) stop("Not all Bands are represented in data file.")
      
      # Identify which rows in ds correspond to the quality data and put into a matrix.
      if(QualityScreen){
        ifelse(any(grepl(QualityBand, ds$row.id)),
               which.QA <- grep(QualityBand, ds$row.id),
               stop("Cannot find quality data in LoadDat. Download quality data with band data in MODISSubsets."))
        QA.time.series <- as.matrix(ds[which.QA,w.ds.dat:ncol(ds)], nrow = length(which.QA), ncol = length(w.ds.dat:ncol(ds)))
      }
      #####
      
      for(bands in 1:length(Bands)){
        
        # Find rows in ds that correspond to Bands[bands] data and put into a matrix.
        ifelse(any(grepl(Bands[bands], ds$row.id)),
               which.band <- grep(Bands[bands], ds$row.id),
               stop("Cannot find band data in LoadDat. Make sure ASCII files in the directory are from MODISSubsets."))
        band.time.series <- as.matrix(ds[which.band,w.ds.dat:ncol(ds)],
                                      nrow = length(which.band), ncol = length(w.ds.dat:ncol(ds)))
        
        ##### Screen pixels in band.time.series: any NoDataFill or pixels < QualityThreshold, will be replaced with NA.
        ifelse(QualityScreen,
               band.time.series <- QualityCheck(Data = band.time.series, QualityScores = QA.time.series, Band = Bands[bands],
                                                NoDataFill = NoDataFill, QualityBand = QualityBand, Product = Product,
                                                QualityThreshold = QualityThreshold),
               band.time.series <- matrix(ifelse(band.time.series != NoDataFill, band.time.series, NA),
                                          nrow = length(which.band)))
        
        # Final check, that band values all fall within the ValidRange (as defined for given MODIS product band).
        if(any(!(band.time.series >= ValidRange[1] && band.time.series <= ValidRange[2]), na.rm = TRUE)){ 
          stop("Some values fall outside the valid range, after no fill values should have been removed.") 
        }
        ####
        
        # Initialise objects for various summaries.
        mean.band <- rep(NA, num.pixels)
        sd.band <- rep(NA, num.pixels)
        band.yield <- rep(NA, num.pixels)
        nofill <- rep(NA, num.pixels)
        poorquality <- rep(NA, num.pixels)
        band.min <- rep(NA, num.pixels)
        band.max <- rep(NA, num.pixels)
       
        if(DiagnosticPlot & num.pixels != 1){
          	cat("Can not provide diagnostic plots for subsets larger than a single pixel")
          	DiagnosticPlot <- FALSE
        }        	
        
        # Run time-series analysis for the ith 
        for(i in 1:ncol(band.time.series)){

          # Minimum and maximum band values observed.
          minobsband <- min(as.numeric(band.time.series[ ,i]) * ScaleFactor, na.rm = TRUE)    
          maxobsband <- max(as.numeric(band.time.series[ ,i]) * ScaleFactor, na.rm = TRUE)
          
          # Assess quality of data at this time-step, to decide how to proceed with analysis.
          data.quality <- ifelse(QualityScreen, sum(!is.na(band.time.series[ ,i])), 2)
          
          if(data.quality >= 2){
            # Linearly interpolate between screened data points, for each pixel, over time (default is daily).
            if(Interpolate){
              if(is.null(InterpolateN)){
                #N <- max(ds$date[!is.na(band.time.series[ ,i])]) - min(ds$date[!is.na(band.time.series[ ,i])])

                InterpolateN <- round(0.9 * as.numeric(length(band.time.series)))
              }
              sout <- approx(x = 1:nrow(band.time.series), y = as.numeric(band.time.series[ ,i]) * ScaleFactor,
                             method = "linear", n = InterpolateN)
            }
            
            # Carry out all the relevant summary analyses, set by options in the function call.
            # (((365*length(years))-16)*365) = average annual yield  (i.e. work out daily yield * 365).
            if(Yield) band.yield[i] <- (sum(sout$y) - minobsband * length(sout$x)) / length(sout$x)
            if(Mean){
              ifelse(Interpolate,
                     mean.band[i] <- mean(sout$y),
                     mean.band[i] <- mean(as.numeric(band.time.series[ ,i]) * ScaleFactor, na.rm = TRUE))
            }
            if(SD){
              ifelse(Interpolate,
                     sd.band[i] <- sd(sout$y),
                     sd.band[i] <- sd(as.numeric(band.time.series[ ,i]) * ScaleFactor, na.rm = TRUE))  
            }
          }
          
          if(data.quality == 1){
            warning("Only single data point that passed the quality screen: cannot summarise", immediate. = TRUE)
            if(Mean) mean.band[i] <- mean(as.numeric(band.time.series[ ,i]) * ScaleFactor, na.rm = TRUE)
          }
          
          if(data.quality == 0) warning("No reliable data for this pixel", immediate. = TRUE)
          
          # Complete final optional summaries, irrespective of data quality.
          if(Min) band.min[i] <- minobsband
          if(Max) band.max[i] <- maxobsband
          
          nofill[i] <- paste(round((sum(ds[ ,i + (w.ds.dat - 1)] == NoDataFill) / length(band.time.series[ ,i])) * 100, 2),
                             "% (", sum(ds[ ,i + (w.ds.dat - 1)] == NoDataFill), "/", length(band.time.series[ ,i]), ")",
                             sep = "")
          
          ifelse(QualityScreen,
                 poorquality[i] <- 
                   paste(round((sum(QA.time.series[ ,i] > QualityThreshold) / length(QA.time.series[ ,i])) * 100, 2),
                         "% (", sum(QA.time.series[ ,i] > QualityThreshold), "/", length(QA.time.series[ ,i]), ")", sep = ""),
                 poorquality[i] <- NA)
       
        } # End of loop for time-series summary analysis for each pixel.
        
        
        # Compile time-series information and relevant summaries data into a final output listed by-sites (.asc files).
        band.data.site[[counter]][(((bands - 1) * num.pixels) + 1):(num.pixels * bands), ] <-
          matrix(data = c(ID = id[ ,counter],
                          lat = rep(lat, num.pixels),
                          long = rep(long, num.pixels),
                          start.date = rep(as.character(min(ds$date)), num.pixels),
                          end.date = rep(as.character(max(ds$date)), num.pixels),
                          data.band = rep(Bands[bands], num.pixels),
                          min.band = band.min,
                          max.band = band.max, 
                          mean.band = mean.band,
                          sd.band = sd.band,
                          band.yield = band.yield,
                          no.fill.data = nofill, 
                          poor.quality.data = poorquality),
                 nrow = num.pixels, ncol = 13)
        
        if(DiagnosticPlot){
        	directory <- file.path(Dir, "DiagnosticPlots")
        	if(file.exists(directory) == FALSE) dir.create(directory)
        	filename <- file.path(directory, paste(id[counter], Product,
                                                 Bands[bands], file.date, ".pdf", sep = "_"))
        	if(data.quality == 1 | data.quality == 0){
          		cat("Too few data points for diagnostic plot for this site\n")
        	} else {
        		pdf(filename)
        		SiteId <- max(paste(lat, long, year))
          		plot(band.time.series[,i] * ScaleFactor, pch = 19, main = id[counter], xlab = "Timestep", ylab = Bands[bands], ylim = c(-0.1, 1))
          		if(Interpolate) lines(sout, col = "red")
          		if(Mean) abline(a = mean.band[i], b=0, lty = 2, col = "red")
          		if(Min) abline(a = band.min[i], b=0, lty = 2)
          		if(Max) abline(a = band.max[i], b=0, lty = 2)
          		dev.off()
        	}
        }
        
      } # End of loop for Bands.
      
    } # End of loop for each ASCII file.
    
    band.data.site <- as.data.frame(do.call("rbind", band.data.site), stringsAsFactors = FALSE)
    names(band.data.site) <- c("ID", "lat", "long", "start.date", "end.date", "data.band", "min.band", "max.band",
                               "mean.band", "sd.band", "band.yield", "no.fill.data", "poor.quality.data")
    
    ### Create SubsetID variable as created in MODISSubsets for correct matching.    
    if(StartDate){
      lat.long <- details[!is.na(details$lat) | !is.na(details$long) | !is.na(details$end.date) | !is.na(details$start.date), ]
      lat.long <- lat.long[!duplicated(data.frame(lat.long$lat, lat.long$long, lat.long$end.date, lat.long$start.date)), ]
    } else if(!StartDate){
      lat.long <- details[!is.na(details$lat) | !is.na(details$long) | !is.na(details$end.date), ]
      lat.long <- lat.long[!duplicated(data.frame(lat.long$lat, lat.long$long, lat.long$end.date)), ]
    }
    
    ID.check <- ifelse(any(names(details) == "ID"), TRUE, FALSE)
    n.unique <- FALSE
    if(ID.check){
      # Check if ID variable can identify unique time series
      n.unique <- length(unique(lat.long$ID)) == nrow(lat.long)
      if(n.unique) names(lat.long)[names(lat.long) == "ID"] <- "SubsetID"
    } else if(!ID.check | !n.unique){
      # Create end.date and start.date
      Year <- FALSE
      POSIXt <- FALSE
      posix.compatible <- try(as.POSIXlt(lat.long$end.date), silent = TRUE)
      if(any(class(lat.long$end.date) == "POSIXt") | all(class(posix.compatible) != "try-error")) POSIXt <- TRUE
      if(all(is.numeric(lat.long$end.date) & nchar(as.character(lat.long$end.date)) == 4) & 
           any(class(posix.compatible) == "try-error")) Year <- TRUE
      if(!Year & !POSIXt) stop("Date information in LoadDat is not recognised as years or as POSIXt format.")
      if(Year & POSIXt) stop("Date information in LoadDat is recognised as both year and POSIXt formats.")
      
      if(Year){
        end.date <- strptime(paste(lat.long$end.date, "-12-31", sep = ""), "%Y-%m-%d")
        if(StartDate){
          start.year.fail <- any(!is.numeric(lat.long$start.date) | nchar(lat.long$start.date) != 4)
          if(start.year.fail) stop("end.date identified as year dates, but start.date does not match.")
          start.date <- strptime(paste(lat.long$start.date, "-01-01", sep = ""), "%Y-%m-%d")
        }
      } else if(POSIXt){
        end.date <- strptime(lat.long$end.date, "%Y-%m-%d")     
        if(StartDate){
          start.posix.fail <- any(class(try(as.POSIXlt(lat.long$end.date), silent = TRUE)) == "try-error")
          if(start.posix.fail) stop("end.date identified as POSIXt dates, but start.date does not match.")
          start.date <- strptime(lat.long$start.date, "%Y-%m-%d")
        }
      }
      if(!StartDate){
        date.ids <- mapply(function(x,y,z) unique(band.data.site$ID[band.data.site$lat == x &
                                                                    band.data.site$long == y &
                                                                    grepl(z, band.data.site$ID)]),
                             x = lat.long$lat, y = lat.long$long, z = as.character(end.date))
        if(!all(sapply(date.ids, length) == 1)){
          stop("Couldn't match summarised data back with original dataset.\n",
               "Make sure LoadDat is exact same dataset used to download MODIS data being summarised.")
        }
        start.date <- substr(date.ids, regexpr("Start", date.ids)+5, regexpr("End", date.ids)-1)
      }
      lat.long <- data.frame(SubsetID = paste("Lat", sprintf("%.5f", lat.long$lat), "Lon", sprintf("%.5f", lat.long$long),
                                              "Start", start.date, "End", end.date, sep = ""),
                             lat.long, stringsAsFactors = FALSE)
    }
    ###
    
    res <- data.frame(cbind(details, matrix(NA, nrow = nrow(details), ncol = 1+(num.pixels * length(Bands)))))
    names(res) <- 
      c(names(details), "SubsetID",
        paste(rep(Bands, each = num.pixels), "_pixel", rep(1:num.pixels, times = length(Bands)), sep = ""))
    
    for(i in 1:length(lat.long$SubsetID)){
      res$SubsetID[sprintf("%.5f", res$lat) %in% sprintf("%.5f", as.numeric(lat.long$lat[i])) &
                     sprintf("%.5f", res$long) %in% sprintf("%.5f", as.numeric(lat.long$long[i])) &
                     res$end.date %in% lat.long$end.date[i]] <- lat.long$SubsetID[i]
      res[res$SubsetID %in% lat.long$SubsetID[i],(which(names(res) == "SubsetID")+1):ncol(res)] <- 
        matrix(band.data.site$mean.band[band.data.site$ID == lat.long$SubsetID[i]],
               nrow = nrow(res[res$SubsetID %in% lat.long$SubsetID[i],(which(names(res) == "SubsetID")+1):ncol(res)]),
               ncol = ncol(res[res$SubsetID %in% lat.long$SubsetID[i],(which(names(res) == "SubsetID")+1):ncol(res)]),
               byrow = TRUE)
    }
    
    #####
    # Write output summary file by appending summary data from all files, producing one file of summary stats output.
    write.table(band.data.site, file = file.path(Dir, paste("MODIS_Summary_", Product, "_", file.date, ".csv", sep = "")),
                sep = ",", row.names = FALSE)
    
    # Write the final appended dataset to a csv file, ready for use, in Dir.
    write.table(res, file = file.path(Dir, paste("MODIS_Data_", Product, "_", file.date, ".csv", sep = "")),
                sep = ",", col.names = TRUE, row.names = FALSE)
    #####
    
    cat("Done! Check the 'MODIS Summary' and 'MODIS Data' output files.\n")
}