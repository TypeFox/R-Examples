MODISTransects <-
function(LoadData, FileSep = NULL, Product, Bands, Size, SaveDir = ".", StartDate = FALSE, TimeSeriesLength = 0)
{
    # Define:  
    # Data are gridded in equal-area tiles in a sinusoidal projection. Each tile consists of a 1200x1200 km data 
    # array of pixels at a finer resolution (see http://modis-land.gsfc.nasa.gov/MODLAND_grid.html).                          
    LONG_EQTR_M = 111.2 * 1000
    
    if(!is.object(LoadData) & !is.character(LoadData)){
      stop("Data is incorrectly specified. Must either be the name of an object in R, or a file path character string.")
    }
    # Load data of locations; external data file, or an R object.
    if(is.object(LoadData)) dat <- data.frame(LoadData)
    if(is.character(LoadData)){
      if(!file.exists(LoadData)) stop("Character string input for LoadData does not resemble an existing file path.")
      if(is.null(FileSep)) stop("Data is a file path, the files delimiter (FileSep) also needed.")
      
      dat <- read.delim(LoadData, sep = FileSep) 
    }
    
    # Check input dataset has variables named as necessary.
    if(!any(names(dat) == "transect") | !any(names(dat) == "start.lat") | 
       !any(names(dat) == "end.lat") | !any(names(dat) == "start.long") | 
       !any(names(dat) == "end.long") | !any(names(dat) == "end.date")){
      stop("Could not find some information that is necessary. May either be missing or incorrectly named.
           See ?MODISTransects for help on data requirements. If data file is loaded, make sure FileSep is sensible.")                                    
    } 
    
    # Check SaveDir matches an existing directory.
    if(!file.exists(SaveDir)) stop("Character string input for SaveDir argument does not resemble an existing file path.")
    
    # Check argument inputs are sensible.
    # If the Product input does not match any product codes in the list output from GetProducts(), stop with error.
    if(!any(Product == GetProducts())) stop("Product entered does not match any available products; see ?GetProducts.")
    
    # If the Bands input does not match with the Product input, stop with error.
    band.test <- lapply(Bands, function(x) !any(x == GetBands(Product)))
    if(any(band.test == TRUE)) stop("A Band does not match Product; see ?GetBands for bands within each product.")
    
    # If Size is not two dimensions or are not integers (greater than expected after rounding, with tolerance around
    # computing precision), stop with error.
    if(!is.numeric(Size)) stop("Size should be numeric class. Two integers.")
    if(length(Size) != 2) stop("Size input must be a vector of integers, with two elements.")
    if(abs(Size[1] - round(Size[1])) > .Machine$double.eps^0.5 | abs(Size[2] - round(Size[2])) > .Machine$double.eps^0.5){
      stop("Size input must be integers.")
    }
    
    if(!is.logical(StartDate)) stop("StartDate must be logical.")
    
    # Year or posixt date format?
    Year <- FALSE
    POSIXt <- FALSE
    
    char.compatible <- as.character(dat$end.date)
    if(!is.character(char.compatible) | all(is.na(char.compatible)) & any(nchar(char.compatible) != 4)) POSIXt <- TRUE
    
    posix.compatible <- try(as.POSIXlt(dat$end.date), silent=TRUE)
    if(class(posix.compatible) == "try-error") Year <- TRUE
    
    if(!Year & !POSIXt) stop("Dates in LoadDat not recognised as years or as POSIXt format.")
    if(Year & POSIXt) stop("Dates in LoadDat recognised as both year and POSIXt formats.")
       
    # Check the start dates are valid.
    if(StartDate){
      if(Year){
        char.compatible <- as.character(dat$start.date)
        if(!is.character(char.compatible) | all(is.na(char.compatible))){
          stop("Year date format detected, but start.date are not compatible with numeric class.")
        }
        if(any(nchar(dat$start.date) != 4)) stop("start.date is not matching year format.")
      } else {
        posix.compatible <- try(as.POSIXlt(dat$start.date), silent = TRUE)
        if(class(posix.compatible) == "try-error") stop("start.date are not all in standard POSIXt format. See ?POSIXt.")
      }    
    }
    
    # Check latitude and longitude inputs are valid coordinate data.
    # Check for missing lat/long data
    if(any(is.na(dat$start.lat) != is.na(dat$start.long) | is.na(dat$start.lat) != is.na(dat$end.lat) | 
           is.na(dat$start.lat) != is.na(dat$end.long) | is.na(dat$start.lat) != is.na(dat$end.date))) { 
      stop("Not equal amount of lats, longs, and dates: there must be locations with incomplete time-series information.") 
    }
    if(abs(dat$start.lat) > 90 || abs(dat$start.long) > 180 || abs(dat$end.lat) > 90 || abs(dat$end.long) > 180){
      stop("Detected some lats or longs beyond the range of valid coordinates.")
    }
    
    if(StartDate){
      if(!any(names(dat) == "start.date")) stop("StartDate == TRUE, but start.date not detected. See ?MODISTransects.")
      if(any(is.na(dat$start.lat) != is.na(dat$start.date))) stop("Not all coordinates have a corresponding start date.")
    } 
    
    # Work out actual width of each pixel in the MODIS projection.
    cell.size.dates <- GetDates(Lat = LoadData$start.lat[1], Long = LoadData$start.long[1], Product = Product)[1:2]
    cell.size <- as.numeric(unname(
      GetSubset(Lat = LoadData$start.lat[1], Long = LoadData$start.long[1], Product = Product, Band = Bands[1],
                StartDate = cell.size.dates[1], EndDate = cell.size.dates[2], KmAboveBelow = 0, KmLeftRight = 0)$pixelsize[[1]]
    ))
    
    # Find all unique transects to download pixels for.
    t.dat <- dat[!duplicated(dat$transect), ]
    cat("Found", nrow(t.dat), "transects. Downloading time-series sets for each transect...\n")
    
    # Loop that reiterates download for each transect.
    for(i in 1:nrow(t.dat)){
      
      # Find the distance, in decimal degrees between the start and end of the transect.
      delta.lat <- t.dat$end.lat[i] - t.dat$start.lat[i]
      delta.long <- round(t.dat$end.long[i] - t.dat$start.long[i], digits = 5)
      # Work out how far in metres is the latitudinal difference between start and end locations.
      lat.metres <- delta.lat * LONG_EQTR_M
      # Determine the curvature angle from the latitude to calculate distance between one longitude at transect location.
      lat.rad <- median(c(t.dat$start.lat[i], t.dat$end.lat[i])) * (pi / 180)
      # Work out how far in metres is the longitudinal difference between start and end locations.
      long.metres <- delta.long * (LONG_EQTR_M * cos(lat.rad))
      # Work out the actual length of the transect.
      transect <- sqrt((lat.metres^2) + (long.metres^2))
      
      # Work out how many points can be equally spaced (i.e. how many pixels) between the start and end coordinates.
      num.points <- transect / cell.size
      # Work out the lat and long distances of the equally spaced points.
      lat.increment <- round(delta.lat / num.points, digits = 5)
      long.increment <- round(delta.long / num.points, digits = 5)
      lat <- t.dat$start.lat[i]
      long <- round(t.dat$start.long[i], digits = 5)
      
      # Take the start lat & long and interpolate new lat & long, using lat & long increments, until the end lat & long.
      # Produces vector of lats and long for all pixels along transect for time-series information in MODISSubsets.
      if(lat.increment > 0){
        
        if(long.increment > 0){                         
          while(lat[length(lat)] <= (t.dat$end.lat[i] - lat.increment) & long[length(long)] <= 
                (round(t.dat$end.long[i], digits=5) - long.increment)){
            lat <- c(lat, round(lat[length(lat)] + lat.increment, digits = 5))
            long <- c(long, round(long[length(long)] + long.increment, digits = 5))
          }
        } else {
          while(lat[length(lat)] <= (t.dat$end.lat[i] - lat.increment) & long[length(long)] >= 
                (round(t.dat$end.long[i], digits=5) - long.increment)){
            lat <- c(lat, round(lat[length(lat)] + lat.increment, digits = 5))
            long <- c(long, round(long[length(long)] + long.increment, digits = 5))
          }
        }
        
      } else {
        
        if(long.increment > 0){
          while(lat[length(lat)] >= (t.dat$end.lat[i] - lat.increment) & long[length(long)] <= 
                (round(t.dat$end.long[i], digits=5) - long.increment)){
            lat <- c(lat, round(lat[length(lat)] + lat.increment, digits = 5))
            long <- c(long, round(long[length(long)] + long.increment, digits = 5))
          }
        } else {
          while(lat[length(lat)] >= (t.dat$end.lat[i] - lat.increment) & long[length(long)] >= 
                  (round(t.dat$end.long[i], digits=5) - long.increment)){
            lat <- c(lat, round(lat[length(lat)] + lat.increment, digits = 5))        
            long <- c(long, round(long[length(long)] + long.increment, digits = 5))   
          }
        }
        
      } # End of if statements that correctly interpolate points along transect line.
      
      # Write vector of end dates & IDs of each transect point to be used for time-series information in MODISSubsets.
      end.date <- rep(t.dat$end.date[i], length(lat))
      ID <- paste("Transect", t.dat$transect[i], "Point", 1:length(lat), sep = "")
      
      # Organise time-series information, with new by-transect IDs for each pixel, for input into MODISSubsets call
      # with optional start date as well as end date.
      if(StartDate){
        start.date <- rep(t.dat$start.date[i], length(lat))                            
        t.subset <- data.frame(ID, lat, long, start.date, end.date)
      } else {
        t.subset <- data.frame(ID, lat, long, end.date)
      }     
      
      #####
      # Do some error checking of pixels in transect before requesting data in MODISSubsets function call.
      xll <- vector(mode = "numeric", length = nrow(t.subset))
      yll <- vector(mode = "numeric", length = nrow(t.subset))
      
      date.for.xy <- GetDates(t.subset$lat[1], t.subset$long[1], Product)[1]
      
      for(n in 1:nrow(t.subset)){
        t.point <- GetSubset(t.subset$lat[n], t.subset$long[n], Product, 
              Bands[1], date.for.xy, date.for.xy, 0, 0)
        xll[n] <- as.numeric(as.character(t.point$xll))
        yll[n] <- as.numeric(as.character(t.point$yll))
      }
      
      # Check which pixels have the same x or y as the previous pixel.
      check.equal.x <- signif(xll[1:length(xll) - 1], digits = 6) == signif(xll[2:length(xll)], digits = 6)
      
      check.equal.y <- signif(yll[1:length(yll) - 1], digits = 5) == signif(yll[2:length(yll)], digits = 5)
      
      # From remaining pixels, check if they are +/- 1 pixel width (i.e. adjacent pixel) away.
      check.new.x <- ifelse(xll[which(!check.equal.x)] < xll[which(!check.equal.x) + 1], 
            round(xll[which(!check.equal.x)]) == round(xll[which(!check.equal.x) + 1] - cell.size), 
            round(xll[which(!check.equal.x)]) == round(xll[which(!check.equal.x) + 1] + cell.size))
      
      check.new.y <- ifelse(yll[which(!check.equal.y)] < yll[which(!check.equal.y) + 1], 
            round(yll[which(!check.equal.y)]) == round(yll[which(!check.equal.y) + 1] - cell.size), 
            round(yll[which(!check.equal.y)]) == round(yll[which(!check.equal.y) + 1] + cell.size))    
      
      if(!all(check.new.x) | !all(check.new.y)){
        # Check if unaccounted for distance between pixels can be attributed to rounding error or MODIS proj uncertainty.
        check.error.x <- ifelse(xll[which(!check.equal.x)] < xll[which(!check.equal.x) + 1], 
              signif(xll[which(!check.equal.x) + 1] - xll[which(!check.equal.x)], digits = 3) == signif(cell.size, digits = 3), 
              signif(xll[which(!check.equal.x) + 1] - xll[which(!check.equal.x)], digits = 3) == -signif(cell.size, digits = 3))
        
        check.error.y <- ifelse(yll[which(!check.equal.y)] < yll[which(!check.equal.y) + 1], 
              signif(yll[which(!check.equal.y) + 1] - yll[which(!check.equal.y)], digits = 3) == signif(cell.size, digits = 3), 
              signif(yll[which(!check.equal.y) + 1] - yll[which(!check.equal.y)], digits = 3) == -signif(cell.size, digits = 3))
        
        # If differences between pixel > expected from checks, abort and produce error stating there are gaps in transect.
        if(any(!check.error.x) | any(!check.error.y)) stop("Error: Gap in transect pixels") 
      }
      
      # Run MODISSubsets to retrieve subset for this transect of pixels.
      MODISSubsets(LoadDat = t.subset, Products = Product, Bands = Bands, Size = Size, SaveDir = SaveDir,
                   StartDate = StartDate, TimeSeriesLength = TimeSeriesLength, Transect = TRUE)
      
    } # End of loop that reiterates download for each transect.
    
  }