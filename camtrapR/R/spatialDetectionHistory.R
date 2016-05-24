spatialDetectionHistory <- function(recordTableIndividual,
                                    species,
                                    camOp,
                                    CTtable,
                                    output,
                                    stationCol = "Station",
                                    speciesCol = "Species",
                                    Xcol,
                                    Ycol,
                                    stationCovariateCols,
                                    individualCol,
                                    individualCovariateCols,
                                    recordDateTimeCol = "DateTimeOriginal",
                                    recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                    occasionLength,
                                    occasionStartTime = 0,
                                    maxNumberDays,
                                    day1,
                                    buffer,
                                    includeEffort = TRUE,
                                    scaleEffort = FALSE,
                                    binaryEffort = FALSE,
                                    timeZone,
                                    makeRMarkInput
)
{

  wd0 <- getwd()
  on.exit(setwd(wd0))

  #################
  # check input
  stopifnot(hasArg(recordTableIndividual))
  stopifnot(hasArg(camOp))

  stopifnot(hasArg(species))
  stopifnot(is.character(species))
  stopifnot(length(species) == 1)

  if(output %in% c("binary", "count") == FALSE) stop("'output' can only be 'binary' or 'count'")

  stopifnot(hasArg(occasionLength))

  stopifnot(c(speciesCol, recordDateTimeCol, stationCol, individualCol) %in% colnames(recordTableIndividual))


  stopifnot(hasArg(stationCol))
  stopifnot(length(stationCol) == 1)
  recordTableIndividual[,stationCol] <- as.character(recordTableIndividual[,stationCol])
  stopifnot(is.character(stationCol))

  stopifnot(hasArg(speciesCol))
  stopifnot(length(speciesCol) == 1)
  recordTableIndividual[,speciesCol] <- as.character(recordTableIndividual[,speciesCol])
  stopifnot(is.character(speciesCol))

  stopifnot(hasArg(individualCol))
  stopifnot(length(individualCol) == 1)
  recordTableIndividual[,individualCol] <- as.character(recordTableIndividual[,individualCol])
  stopifnot(is.character(individualCol))

  stopifnot(hasArg(recordDateTimeCol))
  stopifnot(length(recordDateTimeCol) == 1)
  recordTableIndividual[,recordDateTimeCol] <- as.character(recordTableIndividual[,recordDateTimeCol])
  stopifnot(is.character(recordDateTimeCol))

  if(hasArg(makeRMarkInput)) stopifnot(is.logical(makeRMarkInput))

  stopifnot(c(stationCol, Xcol, Ycol) %in% colnames(CTtable))
  if(any(is.na(CTtable[,Xcol])))stop("there are NAs in Xcol")
  if(any(is.na(CTtable[,Ycol])))stop("there are NAs in Ycol")
  if(any(is.na(CTtable[,stationCol])))stop("there are NAs in stationCol of CTtable")

  if(hasArg(stationCovariateCols)){
    stopifnot(stationCovariateCols %in% colnames(CTtable))
  }

  if(hasArg(timeZone) == FALSE) {
    warning("timeZone is not specified. Assuming UTC", call. = FALSE)
    timeZone <- "UTC"
  }
  if(!is.element(timeZone , OlsonNames())){
    stop("timeZone must be an element of OlsonNames()")
  }


  occasionStartTime    <- as.integer(round(occasionStartTime))
  if(occasionStartTime != 0 & !is.integer(occasionStartTime)) stop ("occasionStartTime must be between 0 and 23")
  if(occasionStartTime < 0 | occasionStartTime >= 24)         stop ("occasionStartTime must be between 0 and 23")

  occasionLength    <- as.integer(round(occasionLength))
  if(occasionLength <= 0)             stop("occasionLength must be a positive integer and not 0")
  if(occasionLength > ncol(camOp)/2)  stop("occasionLength may not be greater than half the total number of days in camOp")
  stopifnot(is.numeric(occasionLength))

  if(hasArg(maxNumberDays)){
    maxNumberDays    <- as.integer(maxNumberDays)
    if(maxNumberDays > ncol(camOp))    stop("maxNumberDays must be smaller than the number of columns of camOp")
    if(maxNumberDays < occasionLength) stop("maxNumberDays must be larger than or equal to occasionLength")
  }

  if(hasArg(buffer)) {
    stopifnot(is.numeric(buffer))
    buffer <- round(buffer)
    stopifnot(buffer >= 1)
  }


  if(species %in% recordTableIndividual[,speciesCol] == FALSE) stop("species is not in speciesCol of recordTableIndividual")

  # check all stations in recordTableIndividual are matched in CTtable
  if(all(recordTableIndividual[,stationCol] %in% CTtable[,stationCol]) == FALSE) {
    stop(paste("items of stationCol in recordTableIndividual are not matched in stationCol of CTtable: ", paste(recordTableIndividual[-which(recordTableIndividual[,stationCol] %in% CTtable[,stationCol]),stationCol], collapse = ", ")))
  }


  if(includeEffort == TRUE){
    if(hasArg(scaleEffort)){
      if(class(scaleEffort) != "logical") stop("scaleEffort must be logical (TRUE or FALSE)")
    } 
    if(hasArg(binaryEffort)){
      if(class(binaryEffort) != "logical") stop("binaryEffort must be logical (TRUE or FALSE)")
    } 
    if(binaryEffort == TRUE & scaleEffort == TRUE) stop("'scaleEffort' and 'binaryEffort' cannot both be TRUE")
  } else {
  scaleEffort  <- FALSE
  binaryEffort <- FALSE
  }


  #####################################################################################################################
  # bring date, time, station ids into shape

  subset_species <- subset(recordTableIndividual, recordTableIndividual[,speciesCol] == species)

  subset_species$DateTime2 <- strptime(subset_species[,recordDateTimeCol], tz = timeZone, format = recordDateTimeFormat)
  subset_species$Date2 <- as.Date(subset_species$DateTime2)

  # check consistency of argument day1
  stopifnot(class(day1) == "character")
  
  if(day1 == "survey") {day1switch <- 1} else {
    if(day1 == "station") {day1switch <- 2} else {
      try(date.test <- as.Date(day1), silent = TRUE)
      if(class(date.test) != "Date") stop("could not interpret argument day1: can only be 'station', 'survey' or a specific date")
      if(all(subset_species$Date2 < as.Date(day1))) stop(paste("there are no records after your specified day1:", day1), call. = FALSE)
      if(hasArg(buffer)) stop("if buffer is defined, day1 can only be 'survey' or 'station'")
      suppressWarnings(rm(date.test))
      day1switch <- 3
    }
  }

  if("POSIXlt" %in% class(subset_species$DateTime2) == FALSE) stop("couldn't interpret recordDateTimeCol of recordTableIndividual using specified recordDateTimeFormat")
  if(any(is.na(subset_species$DateTime2))) stop("at least 1 entry in recordDateTimeCol of recordTableIndividual could not be interpreted using recordDateTimeFormat")

  ####
  cam.op.worked0 <- as.matrix(camOp)

  if(all(as.character(unique(subset_species[,stationCol])) %in% rownames(cam.op.worked0)) == FALSE){
    (stop("Not all values of stationCol in recordTableIndividual are matched by rownames of camOp"))
  }



  ################################################
  # compute date range of stations and records

  cam.tmp.min <- apply(cam.op.worked0, MARGIN = 1, function(X){min(which(!is.na(X)))})    # 1st day of each station
  cam.tmp.max <- apply(cam.op.worked0, MARGIN = 1, function(X){max(which(!is.na(X)))})    # last day of each station

  rec.tmp.min  <- aggregate(subset_species$Date2,
                            list(subset_species[,stationCol]),
                            FUN = min)
  rec.tmp.max  <- aggregate(subset_species$Date2,
                            list(subset_species[,stationCol]),
                            FUN = max)


  date_ranges <- data.frame(rec.min = rec.tmp.min[match(rownames(cam.op.worked0), rec.tmp.min[,1]), 2],
                            rec.max = rec.tmp.max[match(rownames(cam.op.worked0), rec.tmp.max[,1]), 2],
                            cam.min = as.Date(colnames(cam.op.worked0)[cam.tmp.min]),
                            cam.max = as.Date(colnames(cam.op.worked0)[cam.tmp.max])
  )
  date_ranges$diff.days <- date_ranges$cam.max - date_ranges$cam.min
  rownames(date_ranges) <- rownames(cam.op.worked0)

  # check if images were taken between setup and retrieval dates (Error if images outside station date range)
  if(any(date_ranges$rec.min < date_ranges$cam.min | date_ranges$rec.max > date_ranges$cam.max, na.rm = TRUE)) stop(paste("record date outside camera operation date range: ",
                                                                                                                          paste(rownames(date_ranges)[which(date_ranges$rec.min < date_ranges$cam.min)], collapse = ", " )))

  # adjust camera operation dates with buffer and
  # remove records from stations that were operational for a shorter time than the buffer period
  if(hasArg(buffer)){
    if(any(date_ranges$diff.days <= buffer)){
      warning(paste("buffer (", buffer, ") is larger or or equal to date range of cameras: ", paste(rownames(date_ranges)[which(date_ranges$diff.days <= buffer)], collapse = ", ")))
    }
    date_ranges$cam.min <- date_ranges$cam.min + buffer
    if(any(date_ranges$cam.min >= date_ranges$cam.max)){
      stations2remove <- which(date_ranges$cam.min >= date_ranges$cam.max)
      if(length(stations2remove) == nrow(cam.op.worked0)) stop("No station left. your buffer argument is too large.", call. = FALSE)
      if(any(subset_species[,stationCol] %in% rownames(cam.op.worked0) [stations2remove])){
        subset_species <- subset_species[-which(subset_species[,stationCol] %in% rownames(cam.op.worked0) [stations2remove]),]
      }
      if(nrow(subset_species) == 0){stop("No more records left. Choose a smaller buffer argument", call. = FALSE)}
      cam.op.worked0 [stations2remove,] <- NA
      rm(stations2remove)
    }
  }

  rm(cam.tmp.min, cam.tmp.max, rec.tmp.min, rec.tmp.max)


  #######################
  # adjust camera operation matrix

  arg.list0 <- list(cam.op = cam.op.worked0, date_ranges2 = date_ranges, day1_2 = day1, occasionStartTime2 = occasionStartTime)

  if(!hasArg(maxNumberDays) & !hasArg(buffer))  cam.op.worked <- do.call(adjustCameraOperationMatrix, arg.list0)
  if(hasArg(maxNumberDays)  & !hasArg(buffer))  cam.op.worked <- do.call(adjustCameraOperationMatrix, c(arg.list0, maxNumberDays = maxNumberDays))
  if(!hasArg(maxNumberDays) & hasArg(buffer))   cam.op.worked <- do.call(adjustCameraOperationMatrix, c(arg.list0, buffer =  buffer))
  if(hasArg(maxNumberDays)  & hasArg(buffer))   cam.op.worked <- do.call(adjustCameraOperationMatrix, c(arg.list0, maxNumberDays = maxNumberDays, buffer =  buffer))

  rm(arg.list0)


  ######################
  # calculate trapping effort by station and occasion
  effort.tmp <-  calculateTrappingEffort (cam.op = cam.op.worked, occasionLength2 = occasionLength, scaleEffort2 = scaleEffort, includeEffort2 = includeEffort)

  effort <- effort.tmp[[1]]
  if(isTRUE(scaleEffort))     scale.eff.tmp.attr <- effort.tmp[[2]]
  if(isTRUE(binaryEffort))    effort             <- ifelse(effort >= 1, 1, 0)


  ###################
  # remove records that fall into buffer period or were taken after maxNumberDays
  arg.list0 <- list(subset_species2 = subset_species, stationCol2 = stationCol, date_ranges2 = date_ranges, occasionStartTime2 = occasionStartTime, timeZone2 = timeZone)

  if(!hasArg(maxNumberDays) & !hasArg(buffer))  subset_species <- do.call(cleanSubsetSpecies, arg.list0)
  if(hasArg(maxNumberDays)  & !hasArg(buffer))  subset_species <- do.call(cleanSubsetSpecies, c(arg.list0, maxNumberDays = maxNumberDays))
  if(!hasArg(maxNumberDays) & hasArg(buffer))   subset_species <- do.call(cleanSubsetSpecies, c(arg.list0, buffer =  buffer))
  if(hasArg(maxNumberDays)  & hasArg(buffer))   subset_species <- do.call(cleanSubsetSpecies, c(arg.list0, maxNumberDays = maxNumberDays, buffer = buffer))

  rm(arg.list0)


  ############
  #  define the 1st day of the effective survey period.

  if(day1 == "survey")  {    # same starting day for all stations
    time2 = rep(as.POSIXct(paste(as.character(min(date_ranges$cam.min)), "00:00:00"), tz = timeZone), times = nrow(subset_species))
  } else {
    if(day1 == "station") {  # individual starting days for all stations
      time2 = as.POSIXct(paste(as.character(date_ranges$cam.min[match(subset_species[,stationCol], rownames(date_ranges))]), "00:00:00"), tz = timeZone)
    } else {
      if(class(as.Date(day1)) == "Date") {   #  some specific date
        time2 = rep(as.POSIXct( paste(as.Date(day1), "00:00:00"), tz = timeZone), times = nrow(subset_species))
      }
    }
  }


  # calculate time difference between records and first day of detection history (the occasion each record belongs into)
  occasionCol <- "Occasion"

  subset_species[,occasionCol] <- as.numeric(ceiling((difftime(time1 = subset_species[,recordDateTimeCol],
                                                          time2 =  time2,
                                                          units = "secs",
                                                          tz    = timeZone)
                                                 - occasionStartTime * 3600) / (occasionLength * 86400)))



  if(max(subset_species[,occasionCol]) > ncol(effort)) {stop("encountered a bug. Sorry. Please report it.")}


  # column names for effort matrix
  if(includeEffort == TRUE) {colnames(effort) <-  paste("o", seq(1,ncol(effort), by = 1), sep = "")}

  #############
  # build spatial detection history

  # get relevant columns from subset_species
  if(hasArg(individualCovariateCols)){
    sdh0 <- subset_species[, c(individualCol, occasionCol, stationCol, individualCovariateCols)]
    colnames(sdh0) <- c("ID", occasionCol, "trapID", individualCovariateCols)
  } else {
    sdh0 <- subset_species[, c(individualCol, occasionCol, stationCol)]
    colnames(sdh0) <- c("ID", occasionCol, "trapID")
  }

  # remove unidentified individuals
  remove.tmp <- which(is.na(sdh0[,"ID"]))
  if(length(remove.tmp) >= 1){
    sdh0 <- sdh0[-remove.tmp,]
    warning(paste("removed", length(remove.tmp), "records because of missing individual IDs"))
  }

# add required session column
  sdh1 <- data.frame(Session = 1, sdh0)   

# if requested, make remove duplicate records (makes capthist binary. otherwise it returns counts)
  if(output == "binary"){
    sdh2 <- unique(sdh1)                  # remove duplicate rows
  } else {
    sdh2 <- sdh1                          # keep duplicate rows 
  }
  

  ############
  # make secr traps object

  coord.ct <- CTtable[,c(Xcol, Ycol)]
  colnames(coord.ct) <- c("x", "y")
  rownames(coord.ct) <- CTtable[,stationCol]

# set detector type according to desired output (count or binary)
  detectortype <- ifelse (output == "binary", "proximity", "count")

  if(includeEffort == TRUE) {
    if(binaryEffort == TRUE){
      secr.traps <- read.traps(data         = coord.ct,
                               detector     = detectortype,
                               binary.usage = TRUE)
    } else {
      secr.traps <- read.traps(data         = coord.ct,
                               detector     = detectortype,
                               binary.usage = FALSE)
    }
    secr::usage(secr.traps) <- effort
  } else {
    secr.traps <- read.traps(data     = coord.ct,
                             detector = detectortype)
  }

  if(hasArg(stationCovariateCols)){
    whichstationCovariateCols <- which(colnames(CTtable) %in% stationCovariateCols)
    stationCovsDF <- data.frame(CTtable[match(rownames(secr.traps), CTtable[,stationCol]), whichstationCovariateCols])
    if(length(stationCovariateCols) == 1) {
      colnames(stationCovsDF) <- stationCovariateCols
    }
    secr::covariates(secr.traps) <- stationCovsDF
  }

  ###########
  # make capthist object
  capthist.secr <- make.capthist(captures = sdh2,
                                 traps    = secr.traps,
                                 fmt      = "trapID")

if(hasArg(makeRMarkInput)){
  if(isTRUE(makeRMarkInput)){
    RMarkDataframe <- RMarkInput(capthist.secr, covariates = TRUE)
    return(RMarkDataframe)
  }
}

  return(capthist.secr)
}