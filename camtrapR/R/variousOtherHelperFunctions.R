# for functions reading out and tabulating image metadata


addMetadataAsColumns <- function(intable,
                                 metadata.tagname,
                                 metadataHierarchyDelimitor,
                                 multiple_tag_separator)
{
  intable[,metadata.tagname] <- as.character(intable[,metadata.tagname])
  tmp2 <- strsplit(intable[,metadata.tagname], split = ",")      # split items of "HierarchicalSubject" at comma
  tmp3 <- lapply(tmp2, FUN = function(X){X[grep(pattern = metadataHierarchyDelimitor, x = X, fixed = TRUE)]})   # get only the ones with values


  # find all metadata categories, remove spaces
  list.tmp <- vector()
  for(xy in 1:length(tmp3)){
    list.tmp <- c(list.tmp, gsub(pattern = " ",
                                 replacement = "",
                                 x =  unlist(lapply(strsplit(tmp3[[xy]],
                                                             split = metadataHierarchyDelimitor,
                                                             fixed = TRUE),
                                                    FUN = function(Y){Y = Y[1]}))))
  }
  cols2add <- unique(list.tmp)    # these are the columns to add

  # add as columns
  if(length(cols2add) >= 1){    # if anything to add
    intable <- data.frame(intable, matrix(NA, ncol = length(cols2add), nrow = nrow(intable)))
    colnames(intable)[seq((ncol(intable) - length(cols2add) + 1),ncol(intable))] <- cols2add

    # fill metadata columns
    for(xyz in 1:length(cols2add)){
      intable[,cols2add[xyz]] <- unlist(lapply(lapply(tmp3, FUN = function(X) {sapply(strsplit(X[grep(x = X,
                                                                                                      pattern = paste(cols2add[xyz],
                                                                                                                      metadataHierarchyDelimitor,
                                                                                                                      collapse = "",
                                                                                                                      sep      = ""),
                                                                                                      fixed = TRUE)],
                                                                                               split = metadataHierarchyDelimitor,
                                                                                               fixed = TRUE),
                                                                                      FUN = function(Y){Y[2]})}),
                                               FUN = function(Z){paste(Z, collapse = multiple_tag_separator)}))

      intable[which(intable[,cols2add[xyz]] == ""), cols2add[xyz]] <- NA
    } # end for xyz
  } # end if

  which_cols_to_rename <- which(colnames(intable) %in% cols2add)

  # remove spaces and punctuation in column names
  #colnames(intable) <- gsub(pattern = "[[:blank:]]", replacement = "", x = colnames(intable))
  #colnames(intable) <- gsub(pattern = "[[:punct:]]", replacement = "", x = colnames(intable))

  # rename metadata columns with prefix "metadata_"
  colnames(intable)[which_cols_to_rename] <- paste("metadata_", colnames(intable)[which_cols_to_rename], sep = "")

  return(intable)
}



assignSpeciesID <- function(intable,
                            IDfrom,
                            metadataSpeciesTag,
                            speciesCol,
                            dirs_short,
                            i_tmp,
                            multiple_tag_separator)
{

  file.sep <- .Platform$file.sep


  if(IDfrom == "directory"){
    intable[,speciesCol] <-  sapply(strsplit(intable$Directory, split = file.sep, fixed = TRUE), FUN = function(X){X[length(X)]})
    return(intable)
  } else {
    if(hasArg(metadataSpeciesTag)){
      metadataSpeciesTag2 <- paste("metadata", metadataSpeciesTag, sep = "_")
      if(metadataSpeciesTag2 %in% colnames(intable)){

        intable[,speciesCol] <- intable[,metadataSpeciesTag2]
        nrow.intable <- nrow(intable)
        species_records_to_remove <- which(is.na(intable[,speciesCol]))
        if(length(species_records_to_remove) >= 1){
          intable <- intable[-species_records_to_remove,]      #remove records without species tag
          warning(paste( dirs_short[i_tmp],":  removed", length(species_records_to_remove), "records out of", nrow.intable,
                         "because of missing species metadata tag"), call. = FALSE, immediate. = TRUE)
        }

        intable <- separateMultipleSpecies (intable                = intable,
                                            speciesCol             = speciesCol,
                                            multiple_tag_separator = multiple_tag_separator)

        return(intable)
      } else {
        stop(paste("station", dirs_short[i_tmp], ":   metadataSpeciesTag '", metadataSpeciesTag, "' not found in image metadata.", sep = ""), call. = FALSE)
      }
    } else {
      stop(paste("station", dirs_short[i_tmp], ":   cannot figure out species names. Is metadataSpeciesTag defined?"), call. = FALSE)
    }
  }
}


# find and separate multiple species in same image (only if using metadata ID)
separateMultipleSpecies <- function(intable,
                                    speciesCol,
                                    multiple_tag_separator)
{

  records0                 <- intable[,speciesCol]
  records_duplicate        <- strsplit(intable[,speciesCol], split = multiple_tag_separator, fixed = TRUE)
  records_duplicate_length <- sapply(records_duplicate, length)

  if(any(records_duplicate_length > 1)){
    intable <- intable[rep(row.names(intable), records_duplicate_length), ]                # replicate rows with >1 species
    intable[,speciesCol] <- unlist(strsplit (records0, split = multiple_tag_separator))    # assign species anew
  }
  return(intable)
}




# add station and camera id to metadata table

addStationCameraID <- function(intable,
                               dirs_short,
                               stationCol,
                               cameraCol,
                               cameraID,
                               hasStationFolders,
                               i)
{

  file.sep <- .Platform$file.sep

  # append station ID

  if(isTRUE(hasStationFolders)) {       # take station ID from station directories
    intable <- cbind(intable, dirs_short[i])
    colnames(intable)[ncol(intable)] <- stationCol

  } else  {                             # take station ID from image filenames

    station.tmp  <- try(sapply(strsplit(as.character(intable$FileName), split = "__"), FUN = function(X){X[1]}))      # assumes filenames: STATION__Camera__Date/Time(Number).JPG)
    if(length(station.tmp) == nrow(intable)){
      intable <- cbind(intable, station.tmp)
      colnames(intable)[ncol(intable)] <- stationCol
    } else {
      stop(paste(dirs_short[i], ": numbers of images and station ID extracted from image names do not match. Do image filenames begin with station IDs?"))
    }
  }

  # append camera ID

  if(hasArg(cameraID)){
    if(cameraID == "filename"){
      camera.tmp  <- try(sapply(strsplit(as.character(intable$FileName), split = "__"), FUN = function(X){X[2]}))      # assumes filenames: Station__CAMERA__Date/Time(Number).JPG)
      if(length(camera.tmp) == nrow(intable)){
        intable <- cbind(intable, camera.tmp)
        colnames(intable)[ncol(intable)]     <- cameraCol
      }
    }
    if(cameraID == "directory"){            # this can only happen in recordTable. Not in recordTableIndividual
      intable <- cbind(intable,
                       sapply(strsplit(intable$Directory, split = file.sep, fixed = TRUE), FUN = function(X){X[length(X) - 1]}))  # assumes directory structure: Station/Camera/Species
      colnames(intable)[ncol(intable)]     <- cameraCol
    }
  }
  return(intable)
}



#### assess temporal independence between records

assessTemporalIndependence <- function(intable,
                                       deltaTimeComparedTo,
                                       columnOfInterest,
                                       cameraCol,
                                       camerasIndependent,
                                       stationCol,
                                       minDeltaTime)
{

  for(xy in 1:nrow(intable)){     # for every record

    # set independent = TRUE if it is the 1st/only  record of a species / individual
    if(camerasIndependent == TRUE){
      if(intable$DateTimeOriginal[xy]  == min(intable$DateTimeOriginal[which(intable[, columnOfInterest] == intable[xy, columnOfInterest] &
                                                                             intable[, stationCol]       == intable[xy, stationCol] &
                                                                             intable[, cameraCol]        == intable[xy, cameraCol]) ])){    # cameras at same station assessed independently
        intable$independent[xy]       <- TRUE
        intable$delta.time.secs[xy]   <- 0
      }
    } else {
      if(intable$DateTimeOriginal[xy]  == min(intable$DateTimeOriginal[which(intable[, columnOfInterest] == intable[xy, columnOfInterest] &
                                                                             intable[, stationCol]       == intable[xy, stationCol]) ])){
        intable$independent[xy]       <- TRUE
        intable$delta.time.secs[xy]   <- 0
      }
    }

    if(is.na(intable$delta.time.secs[xy])) {   # if not the 1st/only record, calculate time difference to previous records of same species at this station

      if(deltaTimeComparedTo == "lastIndependentRecord"){

        if(camerasIndependent == TRUE){
          which_time2 <- which(intable[, columnOfInterest]       == intable[xy, columnOfInterest] &    # same species/individual
                              intable[, stationCol]              == intable[xy, stationCol] &           # at same station
                              intable[, cameraCol]               == intable[xy, cameraCol] &           # at same camera
                              intable$independent                == TRUE &                             # independent
                              intable$DateTimeOriginal           <  intable$DateTimeOriginal[xy])      # earlier than record xy
        } else {
          which_time2 <- which(intable[, columnOfInterest]       == intable[xy, columnOfInterest] &
                               intable[, stationCol]             == intable[xy, stationCol] &
                               intable$independent               == TRUE &
                               intable$DateTimeOriginal          <  intable$DateTimeOriginal[xy])
        }
      }  else {
        if(camerasIndependent  == TRUE){
          which_time2 <- which(intable[, columnOfInterest]       == intable[xy, columnOfInterest] &
                               intable[, stationCol]             == intable[xy, stationCol] &
                               intable[, cameraCol]              == intable[xy, cameraCol] &
                               intable$DateTimeOriginal          <  intable$DateTimeOriginal[xy])
        } else {
          which_time2 <- which(intable[, columnOfInterest]       == intable[xy, columnOfInterest] &
                               intable[, stationCol]             == intable[xy, stationCol] &
                               intable$DateTimeOriginal          <  intable$DateTimeOriginal[xy])
        }
      }

      #intable$DateTimeOriginal[which_time2] + (minDeltaTime * 60) < intable$DateTimeOriginal[xy]

      # time difference to last (independent) record
      diff_tmp <- min(na.omit(difftime(time1 = intable$DateTimeOriginal[xy],            # delta time to last independent record
                                       time2 = intable$DateTimeOriginal[which_time2],
                                       units = "secs")))

      # save delta time in seconds
      intable$delta.time.secs[xy] <-  diff_tmp
      if(intable$delta.time.secs[xy] >= (minDeltaTime * 60) | intable$delta.time.secs[xy] == 0){
          intable$independent[xy] <- TRUE
        } else {
          intable$independent[xy] <- FALSE
        }

    }   # end   if(intable$DateTimeOriginal[xy] == min(...)} else {...}
  }     # end for(xy in 1:nrow(intable))


  # keep only independent records
  outtable <- intable[intable$delta.time.secs >= (minDeltaTime * 60) |
                      intable$delta.time.secs == 0,]

  return(outtable)
}


# add potential new columns to global record.table

addNewColumnsToGlobalTable <- function(intable,
                                       i,
                                       record.table)
{

  if(i != 1){
    which_cols_to_add_to_d1 <- seq(1, ncol(record.table))[-which(colnames(record.table) %in% colnames(intable))]   # columns in record.table but not in intable

    # if intable lacks columns present in record.table, add them here (filled with NA)
    if(length(which_cols_to_add_to_d1) >= 1){
      intable <- data.frame(intable, as.list(rep(NA, each = length(which_cols_to_add_to_d1))))
      colnames(intable)[(ncol(intable) - length(which_cols_to_add_to_d1) + 1) :  ncol(intable)] <- colnames(record.table)[which_cols_to_add_to_d1]
    }

    # now check which columns are present in intable but not in record.table (new tag groups) and add these (filled with NA)
    which_cols_to_add_to_record.table <- seq(1, ncol(intable))[-which(colnames(intable) %in% colnames(record.table))]  # columns present in intable but not in record.table
    if(length(which_cols_to_add_to_record.table) >= 1){
      record.table <- data.frame(record.table, as.list(rep(NA, each = length(which_cols_to_add_to_record.table))))
      colnames(record.table)[(ncol(record.table) - length(which_cols_to_add_to_record.table) + 1) :  ncol(record.table)] <- colnames(intable)[which_cols_to_add_to_record.table]
    }
    outtable <- intable[,match(colnames(record.table), colnames(intable))]
  } else {
    outtable <- intable
  }
  return(list(outtable, record.table))
}



#####################################################







# for detectionHistory functions


adjustCameraOperationMatrix <- function(cam.op,
                                        date_ranges2,
                                        day1_2,
                                        occasionStartTime2,
                                        maxNumberDays2 = NULL,
                                        buffer2 = NULL){


  #######################
  # adjust camera operation matrix if occasionStartTime2 != 0: make last day NA (it spills into next day in which camera was not set up)

  if(occasionStartTime2 != 0){

    for(k in 1:nrow(cam.op)){
      tmp <- as.vector(cam.op[k,])    # extract camop vector for station k

      tmp[max(which(!is.na(tmp)))]  <- NA                                                     # set last operational day of camOp NA
      tmp[sort(unique(c(grep(pattern = 0, x = tmp), grep(pattern = 0, x = tmp) - 1)))] <- 0   # if station was not operational, set both affected days 0 (the day and the following)
      tmp[which(is.na(as.vector(cam.op[k,])))] <- NA                                          # make sure no additional 0s where there used to be NAs (1st entry)

      cam.op[k,]  <- tmp    # reintroduce that vector into cam op matrix
      rm(tmp)
    }
  }

  ################################################
  #  adjust camera operation matrix if not beginning with 1st day of survey

  # if day1_2 == "survey",  cam.op needs no change because it begins on the 1st day of the survey
  if(day1_2 == "survey"){
    if(!is.null(buffer2)){
      cam.tmp.min <- apply(cam.op, MARGIN = 1, function(X){min(which(!is.na(X)))})

      for(grr in 1:nrow(cam.op)){
        cam.op[grr, seq(1, cam.tmp.min[grr] + buffer2 - 1)] <- NA
      }
      rm(cam.tmp.min)
      cam.op <- cam.op[,-which(apply(cam.op, 2, FUN = function(x){all(is.na(x))}))]   # remove leading NA introduced by buffer2
    }
  } else {
    if(day1_2 == "station") {    # 1st day of each station

      cam.tmp.min <- apply(cam.op, MARGIN = 1, function(X){min(which(!is.na(X)))})
      cam.tmp.max <- apply(cam.op, MARGIN = 1, function(X){max(which(!is.na(X)))})    # last day of each station

      diff.days.tmp <- cam.tmp.max - cam.tmp.min

      cam.op2 <- matrix(NA,
                        nrow = nrow(cam.op),
                        ncol = max(diff.days.tmp)+1)

      # make all stations begin in 1st column
      for(l in 1:nrow(cam.op)){
        cam.op2[l, 1:(diff.days.tmp[l]+1)] <- as.vector(cam.op[l,cam.tmp.min[l]:cam.tmp.max[l]])
      }

      # if buffer2 is defined, remove leading columns
      if(!is.null(buffer2)){
        cam.op2 <- cam.op2[,-seq(1,buffer2)]
      }
      colnames(cam.op2) <- paste("day", 1:ncol(cam.op2), sep = "")
      rownames(cam.op2) <- rownames(cam.op)

      cam.op <- cam.op2
      rm(cam.tmp.min, cam.tmp.max, cam.op2, diff.days.tmp)
    }  else {
      # 1st day = specific date
      if(class(as.Date(day1_2)) == "Date") {
        if(as.Date(day1_2) < min(as.Date(date_ranges2$cam.min)) | as.Date(day1_2) > max(as.Date(date_ranges2$cam.max))) {
          stop(paste("day1_2 (", day1_2, ") is outside the range of days sampled: ", min(as.Date(date_ranges2$cam.min)), " - ", max(as.Date(date_ranges2$cam.max)), sep = ""))
        }
        cam.op <- cam.op [,match(day1_2, colnames(cam.op)) : ncol(cam.op)]
      }
    }
  }


  #######################
  # if maxNumberDays2 is defined, make everything after that NA in cam.op
  if(!is.null(maxNumberDays2)){
    first_day <- apply(cam.op, MARGIN = 1, function(X){min(which(!is.na(X)))})
    for(j in 1:nrow(cam.op)){
      if(first_day[j] + maxNumberDays2 - 1 <= ncol(cam.op)){
        # record.hist   [j, seq(rec.min.tmp[j] + maxNumberDays2, ncol(cam.op), 1)] <- NA
        cam.op [j, seq(from = first_day[j] + maxNumberDays2,
                       to = ncol(cam.op),
                       by = 1)] <- NA
      }
    }
    cam.op <- cam.op[,-which(apply(cam.op, 2, FUN = function(x){all(is.na(x))}))]   # remove tailing NAs introduced by maxNumberDays2
    rm(j, first_day)
  }

  return(cam.op)
}




calculateTrappingEffort <- function(cam.op,
                                    occasionLength2,
                                    scaleEffort2,
                                    includeEffort2){

  ######################
  # calculate trapping effort by station and occasion

  if(occasionLength2 == 1){
    effort <- cam.op          # if occasionLength2 = 1 day, it is identical
  } else {
    effort <- matrix(NA, nrow = nrow(cam.op), ncol = ceiling(ncol(cam.op) / occasionLength2 ))

    index <- 1
    for(m in 1:(ncol(effort))){
      # columns to aggregate
      if(index + occasionLength2 <= ncol(cam.op)){
        index.tmp <- index : (index + occasionLength2 - 1)
      } else {
        index.tmp <- index : ncol(cam.op)
      }
      effort[, m] <- apply(as.matrix(cam.op[,index.tmp]), MARGIN = 1, FUN = sum, na.rm = TRUE)
      effort[, m] <- ifelse(apply(as.matrix(cam.op[,index.tmp]), MARGIN = 1, FUN = function(X) {sum(is.na(X))}) == length(index.tmp), NA, effort[,m])   # if full occasion NA in cam.op, make effort NA
      if(includeEffort2 == FALSE){
        effort[, m] <- ifelse(apply(as.matrix(cam.op[,index.tmp]), MARGIN = 1, FUN = function(X) {any(is.na(X))}), NA, effort[,m])   # if any day NA in cam.op, make effort NA (if effort is not returned)
      }
      index <- index + occasionLength2
    }
    rm(index, index.tmp)
  }


  if(isTRUE(scaleEffort2)){
    if(occasionLength2 == 1) stop("cannot scale effort if occasionLength is 1")
    scale.eff.tmp <- scale(as.vector(effort))
    scale.eff.tmp.attr <- data.frame(effort.scaled.center = NA,
                                     effort.scaled.scale = NA)
    scale.eff.tmp.attr$effort.scaled.center[1] <- attr(scale.eff.tmp, which = "scaled:center")
    scale.eff.tmp.attr$effort.scaled.scale[1]  <- attr(scale.eff.tmp, which = "scaled:scale")
    effort <- matrix(scale.eff.tmp, nrow = nrow(effort), ncol = ncol(effort))
  }

  rownames(effort) <- rownames(cam.op)

  if(isTRUE(scaleEffort2)){
    return(list(effort, scale.eff.tmp.attr))
  } else {
    return(list(effort))
  }
}




cleanSubsetSpecies <- function(subset_species2 ,
                               stationCol2,
                               date_ranges2,
                               buffer2 = NULL,
                               maxNumberDays2 = NULL,
                               occasionStartTime2,
                               timeZone2){

  nrow_subset_species2 <- nrow(subset_species2)

  # remove records that fall into buffer2 period

  if(!is.null(buffer2)){
    corrected_start_time_by_record <- as.POSIXct(paste(as.character(date_ranges2$cam.min[match(subset_species2[,stationCol2], rownames(date_ranges2))]), "00:00:00"), tz = timeZone2) + occasionStartTime2 * 3600
    remove.these <- which(subset_species2$DateTime2 < corrected_start_time_by_record)
    if(length(remove.these) >= 1){
      subset_species2 <- subset_species2[-remove.these,]
      warning(paste(length(remove.these), "records out of", nrow_subset_species2, "were removed because they were before day1"), call. = FALSE, immediate. = TRUE)
      if(nrow(subset_species2) == 0) stop("No more records left. They are all within the buffer period")
      rm(corrected_start_time_by_record, remove.these)
    }
  }

  # remove records that were taken after maxNumberDays2

  if(!is.null(maxNumberDays2)){
    corrected_start_time_by_record <- as.POSIXct(paste(as.character(date_ranges2$cam.min[match(subset_species2[,stationCol2], rownames(date_ranges2))]), "00:00:00"), tz = timeZone2) + occasionStartTime2 * 3600
    remove.these <- which(subset_species2$DateTime2 > corrected_start_time_by_record + maxNumberDays2 * 86400)
    if(length(remove.these) >= 1){
      subset_species2 <- subset_species2[-remove.these,]
      warning(paste(length(remove.these), "records out of", nrow_subset_species2, "were removed because they were after the maxNumberDays"), call. = FALSE, immediate. = TRUE)
      if(nrow(subset_species2) == 0) stop("No more records left.")
      rm(corrected_start_time_by_record, remove.these)
    }
  }

  return(subset_species2)
}







##########################################################################################################











# for function surveyReport

makeSurveyZip <- function(output,
                          recordTable,
                          CTtable ,
                          speciesCol,
                          stationCol,
                          setupCol,
                          retrievalCol,
                          CTDateFormat,
                          CTHasProblems,
                          Xcol,
                          Ycol,
                          recordDateTimeCol,
                          recordDateTimeFormat,
                          sinkpath){

  wd0 <- getwd()
  on.exit(setwd(wd0))

  dir.tmp <- tempdir()
  
  file.sep <- .Platform$file.sep


  dir.zip <- file.path(dir.tmp, paste("surveyReport_", Sys.Date(), sep = ""))
  dir.zip.short <- paste("surveyReport_", Sys.Date(), sep = "")
  unlink(dir.zip, recursive = TRUE)
  dir.create(dir.zip, showWarnings = FALSE, recursive = TRUE)


  # create directories
  invisible(sapply(file.path(dir.zip, c("surveyReport", "activity", "scripts", "detectionMaps")), dir.create, showWarnings = FALSE))


  ######
  # save input tables

  write.csv(recordTable, file = file.path(dir.zip, "recordTable.csv"), row.names = FALSE)
  write.csv(CTtable, file = file.path(dir.zip, "CTtable.csv"), row.names = FALSE)


  ######
  # save surveyReport tables

  dir.tmp2 <- file.path(dir.zip, "surveyReport")

  for(xyz in 1:length(output)){
    write.csv(output[[xyz]],
              file = file.path(dir.tmp2, paste(names(output)[[xyz]], ".csv", sep = "")),
              row.names = FALSE)
  }
  rm(xyz)



  ######
  # make activity plots


  activityDensity(recordTable          = recordTable,
                  allSpecies           = TRUE,
                  speciesCol           = speciesCol,
                  recordDateTimeCol    = recordDateTimeCol,
                  recordDateTimeFormat = recordDateTimeFormat,
                  plotR                = FALSE,
                  writePNG             = TRUE,
                  plotDirectory        = file.path(dir.zip, "activity"))



  ######
  # make detection maps

  if(!is.na(Xcol) & !is.na(Ycol)){
    detectionMaps(CTtable       = CTtable,
                  recordTable   = recordTable,
                  Xcol          = Xcol,
                  Ycol          = Ycol,
                  stationCol    = stationCol,
                  speciesCol    = speciesCol,
                  richnessPlot  = TRUE,
                  speciesPlots  = TRUE,
                  printLabels   = TRUE,
                  plotR         = FALSE,
                  writePNG      = TRUE,
                  plotDirectory = file.path(dir.zip, "detectionMaps"),
                  pngMaxPix     = 1000
    )
  }

  ########################################################################################
  # prepare scripts

  scriptfile <- file.path(dir.zip, "scripts", "camtrapR_scripts.R")
  file.create(scriptfile, showWarnings = FALSE)

  # load basic data

  sink(file = scriptfile)
  cat("###  load data tables  ### \n\n")

  cat("directory.data <- PLEASE_SPECIFY        # this is the directory you got after unzipped the zip file to (e.g. .../surveyReport_2016-02-29/)  \n\n")


  cat("CTtable     <- read.csv(paste(directory.data, 'CTtable.csv', sep = '/'))\n")
  cat("recordTable <- read.csv(paste(directory.data, 'recordTable.csv', sep = '/'))\n\n\n")

  sink()


  # make detection maps   # no, because of coordinate columns

  sink(file = scriptfile, append = TRUE)
  cat("###  plot species detections  ### \n\n")
  if(!is.na(Xcol) & !is.na(Ycol)){
    cat(paste("Xcol <- '", Xcol, "'\n", sep = ""))
    cat(paste("Ycol <- '", Ycol, "'\n\n", sep = ""))
  } else {
    cat("Xcol <- PLEASE_SPECIFY\n")
    cat("Ycol <- PLEASE_SPECIFY\n\n")
  }

  cat(paste("detections <- detectionMaps(CTtable = CTtable,
            recordTable  = recordTable,
            Xcol         = Xcol,
            Ycol         = Ycol,
            stationCol   = '", stationCol, "',
            speciesCol   = '", speciesCol, "',
            writePNG     = FALSE,
            plotR        = TRUE,
            printLabels  = TRUE,
            richnessPlot = TRUE,
            addLegend    = TRUE
  ) \n\n\n", sep = ""))

  sink()


  # camera operation matrix

  sink(file = scriptfile, append = TRUE)
  cat("###  camera operation matrix  ### \n\n")

  cat(paste("cameraOperation <- cameraOperation(CTtable = CTtable,
            stationCol                                  = '", stationCol, "',
            #cameraCol,
            setupCol                                    = '", setupCol, "',
            retrievalCol                                = '", retrievalCol, "',
            hasProblems                                 = '", CTHasProblems, "',
            #byCamera,
            #allCamsOn,
            #camerasIndependent,
            dateFormat                                  = '", CTDateFormat, "' #,
            #writecsv                                   = FALSE,
            #outDir
  ) \n\n\n", sep = ""))

  sink()


  # make detection histories

  sink(file = scriptfile, append = TRUE)
  cat("###  detection histories  ### \n\n")

  cat("day1              <- PLEASE_SPECIFY \n")
  cat("occasionLength    <- PLEASE_SPECIFY\n")
  cat("speciesOfInterest <- PLEASE_SPECIFY\n")
  cat("timeZone          <- PLEASE_SPECIFY \n\n")

  cat(paste("detHist <- detectionHistory(recordTable = recordTable,
            species                                  = speciesOfInterest,
            camOp                                    = cameraOperation,
            stationCol                               = '", stationCol, "',
            speciesCol                               = '", speciesCol, "',
            recordDateTimeCol                        = '", recordDateTimeCol, "',
            recordDateTimeFormat                     = '", recordDateTimeFormat, "',
            occasionLength                           = occasionLength,
            #maxNumberDays,
            day1                                     = day1,
            #buffer                                  = 0,
            includeEffort                            = TRUE,
            scaleEffort                              = FALSE,
            occasionStartTime                        = 0,
            datesAsOccasionNames                     = FALSE,
            timeZone                                 = timeZone,
            writecsv                                 = FALSE #,
            # outDir
  ) \n\n\n", sep = ""))

  sink()


  ### write description text file ###

  readmefile <- file.path(dir.zip, "readme.txt")
  file.create(readmefile, showWarnings = FALSE)

  sink(file = readmefile)


  cat(paste("this zip file contains a summary of a camera trapping survey:\n\n"))
  cat(paste("total survey period:              ", min(output$survey_dates$setup_date), "-", max(output$survey_dates$retrieval_date), "\n"))
  cat(paste("Total number of stations:         ", nrow(output$survey_dates), "\n"))
  cat(paste("Number of operational stations:   ", length(which(output$survey_dates$n_nights_active >= 1)), "\n"))
  cat(paste("Total number of cameras:          ", sum(output$survey_dates$n_cameras), "\n"))
  cat(paste("Total number of active trap days: ", sum(output$survey_dates$n_nights_active), "\n\n"))

  cat("\n-------------------------------------------------------\n\n")

  cat(paste("the following table shows a summary of survey period for each station \n\n"))
  print(output$survey_dates, row.names = FALSE)
  cat("\n-------------------------------------------------------\n\n")

  cat("legend to the data structure in the zip file:\n\n")
  cat(".../activity/              plots of activity density estimations for each species created with activityDensity\n")
  if(!is.na(Xcol) & !is.na(Ycol)){cat(".../detectionMaps/         maps of observed species richness and numbers of species records\n")}
  cat(".../scripts/               a prepared R script for occupancy analyses\n")
  cat(".../surveyReport/          the tables created by surveyReport summarising the survey\n")
  cat(".../CTtable.csv            table of camera trap station IDs, operation times and coordinates\n")
  cat(".../recordTable.csv        table of species records\n")
  cat(".../readme.txt             this file\n\n")

  cat("\n-------------------------------------------------------\n\n")

  cat(paste("species images are located in:\n\n"))
  cat(paste(as.character(recordTable$Directory)[1], "          # the first record\n"))
  cat(paste(as.character(recordTable$Directory)[nrow(recordTable)], "          # the last record\n"))

  cat("\n-------------------------------------------------------\n\n")

  cat("information about who created this report:\n\n")
  print(as.data.frame(Sys.info()))

  cat(paste("\n\n\n\n ***   report created by function surveyReport on",  Sys.time(),  "   ***\n\n\n"))


  sink()

  ######
  # make final zip file

  # list files
  files2zip <- dir(path = dir.zip, full.names = FALSE, recursive = TRUE)
  files2zip <- file.path(dir.zip.short, files2zip)

  # write zip
  setwd(dir.tmp)
  cat("compiling zip file \n",
      paste(sinkpath, paste(dir.zip.short, ".zip\n\n", sep = ""), sep = file.sep))

  suppressMessages(zip(zipfile = file.path(sinkpath,
                                           paste(dir.zip.short, ".zip", sep = "")),
                       files   = files2zip,
                       flags   = ""))



  # remove temporary directory
  unlink(dir.zip, recursive = TRUE)

}