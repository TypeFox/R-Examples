# Function to upload results to a database
# Created 2013 May 15
# Modified 2015 Sept 6

dbUploadAnno <- function(
    annotations,            # Annotations to upload to database: data frame or file path to csv
    survey,                 # Name of survey with which to associate anotations
    db.name='acoustics',    # Connection name in ODBC _and_ on host
    uid,                    # Database User ID, if not in ODBC
    pwd,                    # Database Password, if not in ODBC
    analyst='',             # From `tblPerson`.`pkPersonID`
    ...                     # Additional arguments to RODBC::odbcConnect
){

    if (!requireNamespace("RODBC", quietly = TRUE)) {
        stop("The RODBC package is needed to use this function, but it is not installed. Please install it.", call. = FALSE)
    }  

    start.time <- Sys.time()
    if(any(missing(survey), class(survey) != 'character', length(survey)>1)) stop("Must specify 1 survey name (cannot be a wave object).")
    
    # open the database connection
    if(missing(uid) && missing(pwd)) {dbCon <- RODBC::odbcConnect(db.name, ...)
    } else if(missing(uid)) {dbCon <- RODBC::odbcConnect(db.name, pwd, ...)
    } else dbCon <- RODBC::odbcConnect(db.name, uid, pwd, ...)
    
    # Read in annotations, if necessary
    if(class(annotations) == 'character') { 
      file.ext <- tolower(gsub(".*\\.", "", annotations))
      if(file.ext == 'csv') annotations <- read.csv(annotations) 
      else stop('File extension must be csv, got ', file.ext)
    }     
    # Establish a cleanup procedure
    on.exit(close(dbCon))
    
    # get table of surveys from fldSurveyName
    survey <- RODBC::sqlQuery(dbCon, paste0("Select `pkSurveyID`, `fldOriginalDateModified` FROM `tblSurvey` WHERE `fldSurveyName` = '", survey, "'")) 
        
#    # convert date.time characters to datetime data type format
#    date.time <- unlist(lapply(X=pks.L$date.time, FUN=substr, start=1, stop=19))
#    tzone <- unlist(lapply(pks.L$date.time, function(x) as.character(x, format='%Z')))             

    # the MySQL query to send the hits to the database
    query<- paste0("INSERT INTO `tblAnnotations` (`pkAnnotationID`, `fkSurveyID`, `fkPersonID`, `fldStartTime`, `fldEndTime`, `fldMinFrq`, `fldMaxFrq`, `fldName`) VALUES ('", paste0(NULL, "', '", survey[, 'pkSurveyID'], "', '", analyst, "', '", annotations$start.time, "', '", annotations$end.time, "', '", annotations$min.frq, "', '", annotations$max.frq, "', '", annotations$name, "')", collapse=", ('"))

    # Alert user
    message('\nUploading...')

    # push the query through the open connection to the database     
    status <- RODBC::sqlQuery(dbCon, query)

    # report to user
    message(if(is.na(status[1])) {paste('Done! Upload time:', round(Sys.time()-start.time, 2), 'seconds')
            } else if(status[1] == 'character(0)') {paste('Done! Upload time:', round(Sys.time()-start.time, 2), 'seconds')
            } else paste("Upload unsuccessful; RODBC returned errors: ", paste(status, collapse=" ")))
}    


