# Function to download survey metadata from a database
# Modified 2015 Sept 6

dbDownloadSurvey <- function(
    db.name='acoustics',                # Connection name in ODBC _and_ on host
    uid,                                # Database User ID, if not in ODBC
    pwd,                                # Database Password, if not in ODBC
    start.date,                         # First date to download survey paths 
    end.date,                           # Last date to download survey paths
    loc.prefix,                         # Survey location(s) for which to download survey paths, missing == all
    samp.rate,                          # Specify one or a vector of sampling rate(s)
    ext,                                # Specify one or a vector of extension(s)
    ...                                 # Additional arguments to RODBC::odbcConnect
){

    if (!requireNamespace("RODBC", quietly = TRUE)) {
        stop("The RODBC package is needed to use this function, but it is not installed. Please install it.", call. = FALSE)
    }  
   
    start.time <- Sys.time()
   
    params <- c(!missing(start.date), !missing(end.date), !missing(samp.rate), !missing(ext), !missing(loc.prefix))
    if(sum(params)>0) {
        filter <- " WHERE "
        op <- params
        op[params] <- " AND "
        op[!params] <- ""
    } else {
        filter <- ""
        op <- rep("", 5)
    }
    if(sum(params)<2) {
        filter <- " WHERE "
        op <- rep("", 5)
    }
    
    op <- op[2:5]
    
    if(!missing(loc.prefix)) {
        loc.prefix <- paste0("(", paste0("`tblLocation`.`fldLocationNameAbbreviation` = '", paste0(loc.prefix, collapse="' OR `tblLocation`.`fldLocationNameAbbreviation` = '")), "')")
    } else {
        loc.prefix <- ""
    }
    
    if(!missing(start.date)) {
        start.date <- paste0("(", paste0(" `fldOriginalDateModified` >= '", paste0(start.date, collapse="' OR `tblSurvey`.`fldOriginalDateModified` = '")), "')")
    } else {
        start.date <- ""
    }
    
    if(!missing(end.date)) {
        end.date <- paste0("(", paste0(" `tblSurvey`.`fldOriginalDateModified` <= '", paste0(end.date, collapse="' OR `tblSurvey`.`fldOriginalDateModified` = '")), "')")
    } else {
        end.date <- ""
    }
    
    if(!missing(samp.rate)) {
        samp.rate <- paste0("(", paste0(" `tblSurvey`.`fldSampleRate` = '", paste0(samp.rate, collapse="' OR `tblSurvey`.`fldSampleRate` = '")), "')")
    } else {
        samp.rate <- ""
    }
    
    if(!missing(ext)) {
        ext <- paste0("(", paste0(" `tblSurvey`.`fldRecordingFormat` = '", paste0(ext, collapse="' OR `tblSurvey`.`fldRecordingFormat` = '")), "')")
    } else {
        ext <- ""
    }
    
    query <- paste0("SELECT `fldSurveyName` FROM `tblSurvey` INNER JOIN `tblCardRecorder` ON `tblCardRecorder`.`pkCardRecorderID` = `tblSurvey`.`fkCardRecorderID` INNER JOIN `tblLocation` ON `tblLocation`.`pkLocationID` = `tblCardRecorder`.`fkLocationID`", filter, start.date, op[1], end.date, op[2], samp.rate, op[3], ext, op[4], loc.prefix, ";")
    
    # open the database connection
    if(missing(uid) && missing(pwd)) {dbCon <- RODBC::odbcConnect(db.name, ...)
    } else if(missing(uid)) {dbCon <- RODBC::odbcConnect(db.name, pwd, ...)
    } else dbCon <- RODBC::odbcConnect(db.name, uid, pwd, ...)
    
    # Establish a cleanup procedure
    on.exit(close(dbCon))
    
    # Download the names of the surveys 
    surveys <- RODBC::sqlQuery(dbCon, query, stringsAsFactors=FALSE)
    
    message(if(class(surveys) == 'data.frame') {paste('Done! Download time:', round(Sys.time()-start.time, 2), 'seconds')
            } else paste("Download unsuccessful; RODBC returned errors: ", paste(surveys, collapse=" ")))
    
    return(surveys)
}    



