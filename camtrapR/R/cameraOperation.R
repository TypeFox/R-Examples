cameraOperation <- function(CTtable,
                            stationCol,
                            cameraCol,
                            setupCol,
                            retrievalCol,
                            hasProblems = FALSE,
                            byCamera,
                            allCamsOn,
                            camerasIndependent,
                            dateFormat = "%Y-%m-%d",
                            writecsv = FALSE,
                            outDir){

  # check and prepare input
  wd0 <- getwd()
  on.exit(setwd(wd0))
    
  stopifnot(hasArg(stationCol))
  stopifnot(length(stationCol) == 1)
  CTtable[,stationCol] <- as.character(CTtable[,stationCol])

  stopifnot(hasArg(setupCol))
  stopifnot(length(setupCol) == 1)
  CTtable[,setupCol] <- as.character(CTtable[,setupCol])

  stopifnot(hasArg(retrievalCol))
  stopifnot(length(retrievalCol) == 1)
  CTtable[,retrievalCol] <- as.character(CTtable[,retrievalCol])

  stopifnot(is.logical(writecsv))
  stopifnot(is.logical(hasProblems))

  if(hasArg(byCamera)) {
    stopifnot(is.logical(byCamera))

    if(isTRUE(byCamera) & hasArg(cameraCol) == FALSE){
      stop("if 'byCamera' is TRUE, 'cameraCol' needs to be specified")
    }
  }
  if(hasArg(byCamera) == FALSE & hasArg(cameraCol)){
    stop("if 'cameraCol' is defined, 'byCamera' needs to be specified")
  }


  if(hasArg(cameraCol)){
    if(hasArg(byCamera) == FALSE) stop("if cameraCol is set, byCamera must be specified")
    CTtable[,cameraCol] <- as.character(CTtable[,cameraCol])
    stopifnot(cameraCol %in% colnames(CTtable))
    stopifnot(is.logical(byCamera))
    if(byCamera == FALSE){
      if(hasArg(allCamsOn) == FALSE) stop("if cameraCol is set and byCamera is FALSE, allCamsOn must be specified")
      stopifnot(is.logical(allCamsOn))
      if(allCamsOn == FALSE){
        if(hasArg(camerasIndependent) == FALSE) stop("if cameraCol is set, byCamera is FALSE and allCamsOn is FALSE, camerasIndependent must be specified")
        stopifnot(is.logical(camerasIndependent))
      }
    }
  }


  stopifnot(c(stationCol, setupCol, retrievalCol) %in% colnames(CTtable))

  if(any(is.na(CTtable[,setupCol])))stop("there are NAs in setupCol")
  if(any(is.na(CTtable[,retrievalCol])))stop("there are NAs in retrievalCol")

  if(all(is.na(as.Date(CTtable[,setupCol], format = dateFormat)))){ stop("Cannot read date format in setupCol")}
  if(all(is.na(as.Date(CTtable[,retrievalCol], format = dateFormat)))) {stop("Cannot read date format in retrievalCol")}

  if(any(is.na(as.Date(CTtable[,setupCol], format = dateFormat)))){ stop("at least one entry in setupCol cannot be interpreted using dateFormat")}
  if(any(is.na(as.Date(CTtable[,retrievalCol], format = dateFormat)))) {stop("at least one entry in retrievalCol cannot be interpreted using dateFormat")}

  CTtable[,setupCol] <- as.Date(CTtable[,setupCol], format = dateFormat)
  CTtable[,retrievalCol] <- as.Date(CTtable[,retrievalCol], format = dateFormat)

  if(hasArg(outDir)){
    if(class(outDir) != "character"){stop("outDir must be of class 'character'")}
    if(file.exists(outDir) == FALSE) stop("outDir does not exist")
  }

  if(any(CTtable[,setupCol] > CTtable[,retrievalCol])){
    stop(paste("Setup Date after Retrieval Date:   "),
         paste(CTtable[which(CTtable[,setupCol] > CTtable[,retrievalCol]), stationCol],
               collapse = ", "))
  }

    if(isTRUE(hasProblems)){
         
      cols.prob.from <- grep(colnames(CTtable), pattern = "Problem\\d\\Sfrom")
      cols.prob.to    <- grep(colnames(CTtable), pattern = "Problem\\d\\Sto")
      
      if(length(cols.prob.from) == 0) stop("could not find column ProblemX_from")
      if(length(cols.prob.to) == 0) stop("could not find column ProblemX_to")
          
      if(length(cols.prob.from) != length(cols.prob.to)){
        stop("length of 'Problem..._from' and 'Problem..._to' columns differs. Check format. Sample: 'Problem1_from', 'Problem1_to'")
      }
      
      for(xy in c(cols.prob.from, cols.prob.to)){
        CTtable[,xy] <- as.Date(as.character(CTtable[,xy]), format = dateFormat,  origin = "1970-01-01")
      }

      for(xyz in cols.prob.from){
        if(any(CTtable[,setupCol] > CTtable[,xyz], na.rm = TRUE)){
          stop(paste(paste(CTtable[which(CTtable[,setupCol] > CTtable[,xyz]), stationCol], collapse = ", "), ": Problem begins before Setup"))
        }
      }
      for(xyz2 in cols.prob.to){
        if(any(CTtable[,retrievalCol] < CTtable[,xyz2], na.rm = TRUE)){
          stop(paste(paste(CTtable[which(CTtable[,retrievalCol] < CTtable[,xyz2]), stationCol], collapse = ", "), ": Problem ends after retrieval"))
        }
      }  
      rm(xy, xyz, xyz2)
    }
  
  # function

  if(hasArg(cameraCol)){      # there is a camera column, 1 or more cameras per station
    stationcam <- paste(CTtable[,stationCol], CTtable[,cameraCol], sep = "_")
    m <- matrix(ncol = abs(as.integer(max(CTtable[,retrievalCol]) - min(CTtable[,setupCol]))) + 1,
                nrow  = length(stationcam))
    colnames(m) <- as.character(as.Date(min(CTtable[,setupCol]):max(CTtable[,retrievalCol]), origin = "1970-01-01"))
    rownames(m) <- stationcam

     # fill camera operation matrix (if at least 1 cam operational per Camera)
    unique.tmp <- strsplit(stationcam, split = "_")


    for(i in 1:length(unique.tmp)){

      date0 <- as.character(min(CTtable[,setupCol][CTtable[,stationCol] == unique.tmp[[i]][1] &
                                                   CTtable[,cameraCol]  == unique.tmp[[i]][2]]))

      date1 <- as.character(max(CTtable[,retrievalCol][CTtable[,stationCol] == unique.tmp[[i]][1] &
                                                       CTtable[,cameraCol]  == unique.tmp[[i]][2]]))

      m[i, seq(from = match(date0 , colnames(m)),
               to = match(date1, colnames(m)), by = 1)] <- 1

      if(isTRUE(hasProblems)){
        # set non operational times to 0
        for(j in 1:length(cols.prob.to)){

          date.p0.tmp <- as.character(min(CTtable[CTtable[,stationCol] == unique.tmp[[i]][1] &
                                                  CTtable[,cameraCol]  == unique.tmp[[i]][2], cols.prob.from[j]]))
          date.p1.tmp <- as.character(max(CTtable[CTtable[,stationCol] == unique.tmp[[i]][1] &
                                                  CTtable[,cameraCol]  == unique.tmp[[i]][2], cols.prob.to[j]]))


          if(!is.na(date.p0.tmp) & !is.na(date.p1.tmp)){
            if(date.p1.tmp < date.p0.tmp)stop(paste("Camera", stationcam[i], ", Problem ", j, ": 'to' is smaller than 'from'", sep = ""))
            m[match(stationcam[i], rownames(m)), seq(from = match(date.p0.tmp , colnames(m)),
                                                     to = match(date.p1.tmp, colnames(m)), by = 1)] <- 0
          }
          rm(date.p0.tmp, date.p1.tmp)
        }
      }
      rm(date0, date1)
    }

    if(!isTRUE(byCamera)){  # byCamera = F
      if(isTRUE(allCamsOn)){   # byCamera = F, allCamsOn = T
        dat2 <- aggregate(m, by = list(CTtable[,stationCol]), FUN = min)
        row.names(dat2) <- dat2[,1]
        dat2[,1] <- NULL
        dat2 <- dat2[match(unique(CTtable[,stationCol]), rownames(dat2)),]  # rearrange row order
      } else {     # byCamera = F, allCamsOn = F
        dat2 <- aggregate(m, by = list(CTtable [,stationCol]), FUN = sum, na.rm = TRUE)
        dat2.na <- aggregate(m, by = list(CTtable [,stationCol]), FUN = function(X){all(is.na(X))})
        row.names(dat2) <- row.names(dat2.na) <- dat2[,1]
        dat2[,1] <- dat2.na[,1] <- NULL
        dat2 <- as.matrix(dat2)
        dat2.na <- as.matrix(dat2.na)
        if(any(dat2.na == TRUE)){dat2[which(dat2.na == TRUE)] <- NA}
        rm(dat2.na)

        if(camerasIndependent == FALSE){   
          dat2 <- ifelse(dat2 >= 2,1,dat2) 
        }
        dat2 <- dat2[match(unique(CTtable[,stationCol]), rownames(dat2)),]  # rearrange row order
        dat2 <- as.data.frame(dat2)
      }
    } else {
      dat2 <- as.data.frame(m)
    }
  } else {     # not by camera column, only 1 cam per station

    if(any(table(CTtable[,stationCol]) > 1)){
      stop("at least 1 station has more than 1 item in CTtable. Please specify 'cameraCol'")
    }
    m <- matrix(ncol = abs(as.integer(max(CTtable[,retrievalCol]) - min(CTtable[,setupCol]))) + 1,
                nrow  = length(unique(CTtable[,stationCol])))
    colnames(m) <- as.character(as.Date(min(CTtable[,setupCol]):max(CTtable[,retrievalCol]),
                                        origin = "1970-01-01"))
    rownames(m) <- unique(CTtable[,stationCol])

    unique.tmp <- unique(CTtable[,stationCol])

    for(i in 1:length(unique.tmp)){

      date0 <- as.character(min(CTtable[,setupCol][CTtable[,stationCol] == unique.tmp[i]]))
      date1 <- as.character(max(CTtable[,retrievalCol][CTtable[,stationCol] == unique.tmp[i]]))

      m[match(unique.tmp[i], rownames(m)), seq(from = match(date0 , colnames(m)),
                                               to = match(date1, colnames(m)), by = 1)] <- 1

      # set non operational times to 0
      if(isTRUE(hasProblems)){
        if(length(cols.prob.from) >= 1 & length(cols.prob.to) >= 1){
          for(j in 1:length(cols.prob.to)){
            date.p0.tmp <- as.character(min(CTtable[CTtable[,stationCol] == unique.tmp[i], cols.prob.from[j]]))
            date.p1.tmp <- as.character(max(CTtable[CTtable[,stationCol] == unique.tmp[i], cols.prob.to[j]]))

            if(!is.na(date.p0.tmp) & !is.na(date.p1.tmp)){
              if(date.p1.tmp < date.p0.tmp) stop(paste("Station", unique.tmp[i], ", Problem ", j, ": 'to' is smaller than 'from'", sep = ""))
              if(date.p1.tmp > date1)       stop(paste("Station", unique.tmp[i], ", Problem ", j, ": is outside date range of setup and retrieval", sep = ""))
              if(date.p0.tmp < date0)       stop(paste("Station", unique.tmp[i], ", Problem ", j, ": is outside date range of setup and retrieval", sep = ""))

              m[match(unique.tmp[i], rownames(m)), seq(from = match(date.p0.tmp , colnames(m)),
                                                       to = match(date.p1.tmp, colnames(m)), by = 1)] <- 0
            }
          }
        }
      }
    }
    rm(unique.tmp)
    dat2 <- as.data.frame(m)
  }

  if(writecsv == TRUE){

    hasProblemsString <- ifelse(isTRUE(hasProblems), "with_problems_", "")
    byCameraString <- ifelse(isTRUE(byCamera), "by_camera", "by_station")
    filename.out <- paste("CameraOperationMatrix_",byCameraString, "_", hasProblemsString, Sys.Date(), ".csv", sep = "")

    if(hasArg(outDir) == FALSE){
      setwd(getwd())
      write.csv(dat2, file = filename.out,
                row.names = TRUE)
    } else {
      setwd(outDir)
      write.csv(dat2, file = filename.out,
                row.names = TRUE)
    }
  }
  return(as.matrix(dat2))
}