detectionHistory <- function(recordTable,
                             species,
                             camOp,
                             stationCol = "Station",
                             speciesCol = "Species",
                             recordDateTimeCol = "DateTimeOriginal",
                             recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                             occasionLength,
                             maxNumberDays,
                             day1,
                             buffer,
                             includeEffort = TRUE,
                             scaleEffort = FALSE,
                             occasionStartTime = 0,
                             datesAsOccasionNames = FALSE,
                             timeZone,
                             writecsv = FALSE,
                             outDir)
{
  wd0 <- getwd()
  on.exit(setwd(wd0))
  #################
  # check input

  stopifnot(hasArg(species))
  stopifnot(is.character(species))
  stopifnot(length(species) == 1)

  stopifnot(hasArg(occasionLength))

  stopifnot(hasArg(recordTable))
  stopifnot(hasArg(camOp))

  stopifnot(hasArg(stationCol))
  stopifnot(length(stationCol) == 1)
  recordTable[,stationCol] <- as.character(recordTable[,stationCol])
  stopifnot(is.character(stationCol))


  stopifnot(hasArg(speciesCol))
  stopifnot(length(speciesCol) == 1)
  recordTable[,speciesCol] <- as.character(recordTable[,speciesCol])
  stopifnot(is.character(speciesCol))


  stopifnot(hasArg(recordDateTimeCol))
  stopifnot(length(recordDateTimeCol) == 1)
  recordTable[,recordDateTimeCol] <- as.character(recordTable[,recordDateTimeCol])
  stopifnot(is.character(recordDateTimeCol))


  if(hasArg(timeZone) == FALSE) {
    warning("timeZone is not specified. Assuming UTC", call. = FALSE)
    timeZone <- "UTC"
  }
  if(!is.element(timeZone , OlsonNames())){
    stop("timeZone must be an element of OlsonNames()")
  }
  stopifnot(is.logical(writecsv))

  occasionStartTime <- as.integer(round(occasionStartTime))
  if(occasionStartTime != 0 & !is.integer(occasionStartTime)) {stop ("occasionStartTime must be between 0 and 23")}
  if(occasionStartTime < 0 | occasionStartTime >= 24){stop ("occasionStartTime must be between 0 and 23")}

  occasionLength <- as.integer(round(occasionLength))
  if(occasionLength <= 0) stop("occasionLength must be a positive integer and not 0")
  if(occasionLength > ncol(camOp)/2) stop("occasionLength may not be greater than half the total number of days in camOp")
  stopifnot(is.numeric(occasionLength))

  if(hasArg(maxNumberDays)){
    maxNumberDays <- as.integer(maxNumberDays)
    if(maxNumberDays > ncol(camOp)) stop("maxNumberDays must be smaller than the number of columns of camOp")
    if(maxNumberDays < occasionLength) stop("maxNumberDays must be larger than or equal to occasionLength")
  }
  
  if(hasArg(buffer)) {
    stopifnot(is.numeric(buffer))
    buffer <- round(buffer)
    stopifnot(buffer >= 1)
  }


  stopifnot(c(speciesCol, recordDateTimeCol, stationCol) %in% colnames(recordTable))


  if(species %in% recordTable[,speciesCol] == FALSE) stop("species is not in speciesCol of recordTable")

  if(writecsv == TRUE){
    if(file.exists(outDir) == FALSE){stop("outDir does not exist")}
  }


  if(includeEffort == TRUE){
    if(!hasArg(scaleEffort)) stop("scaleEffort must be defined if includeEffort is TRUE")
    if(class(scaleEffort) != "logical") stop("scaleEffort must be logical (TRUE or FALSE)")
  } else {scaleEffort <- FALSE}


  #############
  # bring date, time, station ids into shape

  subset_species           <- subset(recordTable, recordTable[,speciesCol] == species)
  subset_species$DateTime2 <- strptime(subset_species[,recordDateTimeCol], tz = timeZone, format = recordDateTimeFormat)
  subset_species$Date2     <- as.Date(subset_species$DateTime2)


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

  if("POSIXlt" %in% class(subset_species$DateTime2) == FALSE) stop("couldn't interpret recordDateTimeCol of recordTable using specified recordDateTimeFormat")
  if(any(is.na(subset_species$DateTime2))) stop("at least 1 entry in recordDateTimeCol of recordTable could not be interpreted using recordDateTimeFormat")


  ####
  cam.op.worked0 <- as.matrix(camOp)

  if(all(as.character(unique(subset_species[,stationCol])) %in% rownames(cam.op.worked0)) == FALSE){
    (stop("Not all values of stationCol in recordTable are matched by rownames of camOp"))
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


  date_ranges <- data.frame(rec.min = rec.tmp.min[match(rownames(cam.op.worked0), rec.tmp.min[,1]), 2],   # as.Date(unlist(rec.tmp.min), origin = "1970-01-01"),
                            rec.max = rec.tmp.max[match(rownames(cam.op.worked0), rec.tmp.max[,1]), 2],   # as.Date(unlist(rec.tmp.max), origin = "1970-01-01"),
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
  if(isTRUE(scaleEffort))  scale.eff.tmp.attr <- effort.tmp[[2]]


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

  subset_species$occasion <- as.numeric(ceiling((difftime(time1  = subset_species[,recordDateTimeCol],
                                                          time2 =  time2,
                                                          units = "secs",
                                                          tz = timeZone)
                                                 - occasionStartTime * 3600) / (occasionLength * 86400)))


  ############
  # make detection history

  if(occasionLength == 1){
    record.hist <- cam.op.worked
    record.hist <- ifelse(record.hist == 0, NA, record.hist)    # if cameras not operational, set NA
    record.hist <- ifelse(record.hist >= 1, 0,  record.hist)    # if cameras operational, set to 0

    occasions.by.station <- tapply(X = subset_species$occasion, INDEX = subset_species[,stationCol], FUN = unique, simplify = FALSE)

    for(xyz in which(sapply(occasions.by.station, FUN = function(x){!is.null(x)}))){
      record.hist[match(names(occasions.by.station)[xyz], rownames(record.hist)), occasions.by.station[[xyz]]] <- 1
    }
    record.hist[is.na(cam.op.worked)] <- NA   # remove the records that were taken when cams were NA (redundant with above:   # remove records taken after day1 + maxNumberDays)

    rm(occasions.by.station, xyz)
  } else {                                    # occasion lenght > 1
    record.hist <- effort
    record.hist <- ifelse(!is.na(record.hist), 0, record.hist)
    rownames(record.hist) <- rownames(cam.op.worked)

    occasions.by.station <- tapply(X = subset_species$occasion, INDEX = subset_species[,stationCol], FUN = unique, simplify = FALSE)

    for(xyz in which(sapply(occasions.by.station, FUN = function(x){!is.null(x)}))){
      record.hist[match(names(occasions.by.station)[xyz], rownames(record.hist)), occasions.by.station[[xyz]]] <- 1
    }
    record.hist[is.na(effort)] <- NA
  }



  # rownames of output
  row.names(record.hist) <- row.names(effort) <- rownames(cam.op.worked)

  # column names for output table
  colnames(cam.op.worked) <- paste(colnames(cam.op.worked), "+", occasionStartTime, "h", sep = "")

# assign column names
  if(isTRUE(datesAsOccasionNames)){
    if(occasionLength == 1){
      colnames.tmp <- colnames(cam.op.worked)
    } else {
      seq.tmp <- seq(from = 1, by = occasionLength, length.out = ncol(record.hist))
      colnames.tmp <- paste(colnames(cam.op.worked)[seq.tmp],
                            colnames(cam.op.worked)[seq.tmp + occasionLength - 1], sep = "_")
      colnames.tmp[length(colnames.tmp)] <- paste(colnames(cam.op.worked)[max(seq.tmp)],      # fix the last one
                                                  colnames(cam.op.worked)[ncol(cam.op.worked)],
                                                  sep = "_")
    }
    colnames(record.hist) <- colnames(effort) <- colnames.tmp

  }  else {
    colnames(record.hist) <- colnames(effort) <-  paste("o", seq(1,ncol(record.hist), by = 1), sep = "")
  }

  ################################################
  # save output as table
  if(day1switch == 1) day1string <- "_first_day_from_survey"
  if(day1switch == 2) day1string <- "_first_day_by_station"
  if(day1switch == 3) day1string <- paste("_first_day", day1, sep = "_")

  effortstring <- ifelse(isTRUE(includeEffort), "with_effort__", "no_effort__")
  maxNumberDaysstring <- ifelse(hasArg(maxNumberDays), paste("max",maxNumberDays,"days_", sep = ""), "")
  if(isTRUE(includeEffort)){
    scaleEffortstring <- ifelse(isTRUE(scaleEffort), "scaled_", "not_scaled_")
  } else {
    scaleEffortstring <- ""
  }

  outtable.name <- paste(species, "__record_history__", effortstring,
                         occasionLength, "_days_per_occasion_",
                         maxNumberDaysstring,
                         "_occasionStart", occasionStartTime,"h_",
                         day1string, "__",
                         Sys.Date(),
                         ".csv", sep = "")

  outtable.name.effort <- paste(species, "__effort__",
                                scaleEffortstring,
                                occasionLength, "_days_per_occasion_",
                                maxNumberDaysstring,
                                "_occasionStart", occasionStartTime,"h_",
                                day1string, "__",
                                Sys.Date(),
                                ".csv", sep = "")

  outtable.name.effort.scale <- paste(species, "__effort_scaling_parameters__",
                                      occasionLength, "_days_per_occasion_",
                                      maxNumberDaysstring,
                                      "_occasionStart", occasionStartTime,"h_",
                                      day1string, "__",
                                      Sys.Date(),
                                      ".csv", sep = "")

  if(isTRUE(writecsv)){
    setwd(outDir)
    write.csv(record.hist, file = outtable.name)
    if(isTRUE(includeEffort)){
      write.csv(effort, file = outtable.name.effort)
      if(hasArg(scaleEffort)){
        if(scaleEffort == TRUE)  write.csv(scale.eff.tmp.attr, file = outtable.name.effort.scale)
      }
    }
  }
  if(isTRUE(includeEffort)){
    if(scaleEffort == TRUE){
      return(list(detection_history = record.hist,
                  effort = effort,
                  effort_scaling_parameters = scale.eff.tmp.attr))
    } else {
      return(list(detection_history = record.hist,
                  effort = effort))
    }
  } else {
    return(list(detection_history = record.hist))
  }
}