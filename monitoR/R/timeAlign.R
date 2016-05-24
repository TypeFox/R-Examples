# Function to condense detections or peaks from a number of templates (of the same species and detection type); events that occur within an adjustable time buffer of one another are considered to be duplicated.  In such cases only the event that had the highest score is saved.  The buffer is only useful when the signal is centered in the template!  (could include a centering function for accuracy)  Input is an object of class'detectionList', a single data frame, or list of either file paths or data frames.  Returns a single data frame.  Must be done separately for each survey.

# Written 2013 Mar 29
# Modified 2014 MAR 29

timeAlign <- function(
  x,                 # Object of class'detectionList', single data frame with results from multiple templates, list of file paths or data frames with results from one or more templates.
  what='detections', # If an detection list is provided, should detections or peaks be used?
  tol=1              # Alignment checks to see if the surrounding events are within 0.5 * this value, in seconds 
  ) {  

    # cut the tolerance value in half to produce a +/- time buffer
    tol <- 0.5*tol
    # character class items are assumed to be file paths to csv files...check
    if(class(x) == 'detectionList' && tolower(what) %in% c('detections', 'd', 'det')) {
      cat('Getting detections\n')
      all.results <- getDetections(x)
    } else if(class(x) == 'detectionList' && tolower(what) %in% c('peaks', 'p', 'pks')) {
      cat('Getting peaks\n')
      all.results <- getPeaks(x)
    } else if(class(x) == 'data.frame') {all.results <- x
    } else if(class(x) == 'list') {
      check.class <- unlist(lapply(as.list(x), function(file) class(file)))
        if(all(check.class == 'character')) {
            all.results <- lapply(x, function(file) {
              chars <- nchar(file)
              ext <- tolower(gsub(".*\\.", "", file))
                if(all(ext == 'csv')) {all.results <- lapply(x, function(data) read.csv(data, stringsAsFactors=FALSE))
                } else stop('x must be a csv file.')
                return(all.results)
                })
        } else if(all(check.class == 'data.frame')) all.results <- rbindf(x)
    } else stop('x is unfamiliar; must be detectionList, data frame, or a list of either data frames or file paths to csv files.')
    # re-order data frames to prepare to merge them
    colNm <- names(all.results)
    all.results <- all.results[, colNm]
    # Determine source of x; must be from findPeaks or from acoustics database
    names.res <- names(all.results)
     if(all(c("fldTime", "fldTemplateName", "fldScore") %in% names.res)) {
        time.fld <- "fldTime"
        temp.fld <- "fldTemplateName"
        scor.fld <- "fldScore"
        onAmp.fld <- "fldOnAmp"
        offAmp.fld <- "fldOffAmp"
    } else if(all(c("time", "template", "score") %in% names.res)) {
        time.fld <- "time"
        temp.fld <- "template"
        scor.fld <- "score"
        onAmp.fld <- "on.amp"
        offAmp.fld <- "off.amp"
    } else stop("Unrecognized column names.  x must be from findPeaks or from an acoustics database.")
    # Determine whether to include the amp stuff
    if(any(c("fldOnAmp", "fldOffAmp", "on.amp", "off.amp") %in% names.res)) {report.amp <- TRUE
    } else report.amp <- FALSE
    # Order by time
    all.results <- all.results[order(all.results[, time.fld]), ]
    row.names(all.results) <- 1:nrow(all.results)
    # Determine if the time in each row is the same (+/- 1 second) as the row above
#    cat('Identifying duplicates\n')
    nRows <- nrow(all.results)
    earlytimes <- all.results[c(1:(nRows-1)), time.fld]
    latetimes <- all.results[c(2:nRows), time.fld]
    same1below <- c(earlytimes+tol >= latetimes & earlytimes+tol <= latetimes+tol, FALSE)
    same1above <- c(FALSE, earlytimes+tol >= latetimes & earlytimes+tol <= latetimes+tol)
    group <- same1above*same1below
    # prepare to align fields
    nTimes <- nScore <- nOnAmp <- nOffAmp <- nTemplateName <- rep(NA, nRows)
    # identify groups of 3 or more, assign them all data from line with the max score
    groupID <- rep(NA, nRows)
    gid <- 1
    for(i in 2:(nRows-1)) {
      if(group[i]) {groupID[(i-1):(i+1)] <- gid
      } else if(!group[i] & group[i-1]) gid <- gid+1
    }
    
    for(i in na.omit(unique(groupID))) {
      rID <- which(groupID == i)
      num.evt <- length(rID)
      rmaxscor <- rID[which(all.results[rID, scor.fld] == max(all.results[rID, scor.fld]))]
      if(length(rmaxscor)>1) rmaxscor <- rmaxscor[1]
      nScore[rID] <- rep(all.results[rmaxscor, scor.fld], num.evt)
      nTimes[rID] <- rep(all.results[rmaxscor, time.fld], num.evt)
      if(report.amp) {
          nOnAmp[rID] <- rep(all.results[rmaxscor, onAmp.fld], num.evt)  
          nOffAmp[rID] <- rep(all.results[rmaxscor, offAmp.fld], num.evt)
      }
      nTemplateName[rID] <- rep(all.results[rmaxscor, temp.fld], num.evt)
    }
    # identify groups of 2, assign one data from line with the max score
    groupID <- rep(NA, nRows)
    gid <- 0
    for(i in 1:nRows) {
      if(any(same1above[i], same1below[i]) && is.na(nScore[i])) {
          if(same1above[i]) groupID[i] <- gid
          else {
              gid <- gid+1
              groupID[i] <- gid
          }
      } #else if(all(!same1above[i], !same1below[i], !group[i])) gid <- gid+1
    }
    for(i in na.omit(unique(groupID))) {
      rID <- which(groupID == i)
      num.evt <- length(rID)
      rmaxscor <- rID[which(all.results[rID, scor.fld] == max(all.results[rID, scor.fld]))]
      if(length(rmaxscor)>1) rmaxscor <- rmaxscor[1]
      nScore[rID] <- rep(all.results[rmaxscor, scor.fld], num.evt)
      nTimes[rID] <- rep(all.results[rmaxscor, time.fld], num.evt)
      if(report.amp) {
          nOnAmp[rID] <- rep(all.results[rmaxscor, onAmp.fld], num.evt)
          nOffAmp[rID] <- rep(all.results[rmaxscor, offAmp.fld], num.evt)
      }
      nTemplateName[rID] <- rep(all.results[rmaxscor, temp.fld], num.evt)
    }
    # leave remaining events unchanged
    rID <- which(is.na(nScore))
    nScore[rID] <- all.results[rID, scor.fld]
    nTimes[rID] <- all.results[rID, time.fld]
    if(report.amp) {
        nOnAmp[rID] <- all.results[rID, onAmp.fld]
        nOffAmp[rID] <- all.results[rID, offAmp.fld]
    }
    nTemplateName[rID] <- all.results[rID, temp.fld]
    # Reassemble data with new fields
    all.results[time.fld] <- nTimes
    all.results[scor.fld] <- nScore
    if(report.amp) {
        all.results[onAmp.fld] <- nOnAmp
        all.results[offAmp.fld] <- nOffAmp
    }
    all.results[temp.fld] <- nTemplateName
    # filter for unique values
#    cat('Dropping lower scoring duplicates\n')
    goodRows <- cbind(all.results[time.fld], all.results[scor.fld])
    goodRows <- unique(goodRows)
    goodRows <- as.numeric(row.names(goodRows))
    all.results <- all.results[goodRows, ]
    row.names(all.results) <- 1:nrow(all.results)
#    #save the original template names
#    all.results['fldOrigTemplate'] <- all.results[temp.fld]
#    all.results[temp.fld] <- 'all'
#    # Add all columns to the main data frame
#    all.results['manual'] <- manual
#    all.results['sameabove'] <- same1above
#    all.results['samebelow'] <- same1below
#    # Return the data
    return(all.results)
  }
