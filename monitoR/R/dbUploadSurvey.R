# Function to upload survey metadata to a database, also updates existing metadata
# Created 2012 Oct 23
# Modified 2015 Sept 6

dbUploadSurvey <- function(
    db.name='acoustics',            # Connection name in ODBC _and_ on host
    uid,                            # Database User ID, if not in ODBC
    pwd,                            # Database Password, if not in ODBC
    survey.meta,                    # Survey metadata to upload to database
    update.query=FALSE,             # FALSE -> INSERT INTO query, TRUE -> UPDATE-WHERE query
    tz,                             # Time zone, if not specified in file name or metadata
    ...                             # Additional arguments to RODBC::odbcConnect
    ){

    if (!requireNamespace("RODBC", quietly = TRUE)) {
        stop("The RODBC package is needed to use this function, but it is not installed. Please install it.", call. = FALSE)
    }  
    
    start.time <- Sys.time()
  
    # open the database connection
    if(missing(uid) && missing(pwd)) {dbCon <- RODBC::odbcConnect(db.name, ...)
    } else if(missing(uid)) {dbCon <- RODBC::odbcConnect(db.name, pwd, ...)
    } else dbCon <- RODBC::odbcConnect(db.name, uid, pwd, ...)
    # Establish a cleanup procedure
    on.exit(close(dbCon))
    
    # Unpack modification times
    if(!'fldOriginalDateModified' %in% names(survey.meta)) {
        date.time.info <- regmatches(survey.meta[, 'fldSurveyName'], regexpr('[0-9]{4}-[0-9]{2}-[0-9]{2}[ _][0-9]{6}[ _][A-Z][SDM]T', survey.meta[, 'fldSurveyName']))
        if(length(date.time.info)>=1 && nchar(date.time.info) == 21) {
            dm <- substr(survey.meta[, 'fldSurveyName'], start=8, stop=24)
            dm <- as.Date(dm, format="%Y-%m-%d_%H%M%S")
            dm <- as.POSIXct(dm, tz=tz, format="%Y-%m-%d_%H%M%S %Z")
            survey.meta['fldOriginalDateModified'] <- as.character(dm, format="%Y-%m-%d_%H%M%S %Z")
        } else stop('No "fldOriginalDateModified" in metadata and file name does not have date and time info.')
        }
	if(!update.query) { 
		mtimes <- survey.meta[, 'fldOriginalDateModified']
		dates <- substr(mtimes, start=1, stop=11)
		hh <- substr(mtimes, start=12, stop=13)
		mm <- substr(mtimes, start=14, stop=15)
		ss <- substr(mtimes, start=16, stop=17)
		date.time <- paste0(unlist(dates), unlist(hh), ":", unlist(mm), ":", unlist(ss))
		tzone <- substr(mtimes, start=19, stop=21)
		}
	# Fill out missing fields
	empty <- !c('fkCardRecorderID', 'fldSurveyLength', 'fldOriginalDateModified', 'fldTimeZone', 'fldOriginalRecordingName', 'fldSurveyName', 'fldRecordingFormat', 'fldSampleRate', 'fldBitsperSample', 'fldChannels') %in% names(survey.meta) 
	empty <- c('fkCardRecorderID', 'fldSurveyLength', 'fldOriginalDateModified', 'fldTimeZone', 'fldOriginalRecordingName', 'fldSurveyName', 'fldRecordingFormat', 'fldSampleRate', 'fldBitsperSample', 'fldChannels')[empty]
	survey.meta[empty] <- NA
    # test conditions and write out the possible MySQL queries
    if(!update.query) { 
        # the MySQL query to send the survey metadata to the database
        query <- paste("INSERT INTO `tblSurvey` (`pkSurveyID`, `fkCardRecorderID`, `fldSurveyLength`, `fldOriginalDateModified`, `fldTimeZone`, `fldOriginalRecordingName`, `fldSurveyName`, `fldRecordingFormat`, `fldSampleRate`, `fldBitsperSample`, `fldChannels`) VALUES ('", paste(NULL, "', '", survey.meta$fkCardRecorder, "', '", survey.meta$fldSurveyLength, "', '", date.time, "', '", tzone, "', '", survey.meta$fldOriginalRecordingName, "', '", survey.meta$fldSurveyName, "', '", survey.meta$fldRecordingFormat, "', '", survey.meta$fldSampleRate, "', 	'", survey.meta$fldBitsperSample, "', '", survey.meta$fldChannels, "')", sep="", collapse=", ('"), sep="")
    
    } else {
        # the MySQL query to send to update the survey metadata        
        query <- paste("UPDATE `tblSurvey` SET `fldSurveyLength` = '", survey.meta$fldSurveyLength, "', `fldSampleRate` = '", survey.meta$fldSampleRate, "', `fldBitsperSample` = '", survey.meta$fldBitsperSample, "', `fldChannels` = '", survey.meta$fldChannels, "' WHERE `fldSurveyName` = '", survey.meta$fldSurveyName, "'", sep="")
    }
    
    # push the query through the open connection to the database    
    if(update.query) {status <- lapply(X=query, FUN=RODBC::sqlQuery, channel=dbCon)
    } else {status <- RODBC::sqlQuery(dbCon, query)}
    
    # report to user
    message(if(is.na(status[1])) {paste('Done! Upload time:', round(Sys.time()-start.time, 2), 'seconds')
            } else if(status[1] == 'character(0)') {paste('Done! Upload time:', round(Sys.time()-start.time, 2), 'seconds')
            } else paste("Upload unsuccessful; RODBC returned errors: ", paste(status, collapse=" ")))
}



