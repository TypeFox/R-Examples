# Function to download card and recorder id by location and date from database
# Modified 2015 Sept 6

dbDownloadCardRecorderID <- function(
    db.name='acoustics',        # Connection name in ODBC _and_ on host
    uid,                        # Database User ID, if not in ODBC
    pwd,                        # Database Password, if not in ODBC
    date.deployed,              # First deploy date to download 
    date.collected,             # Last collection date to download
    loc.prefix,                 # Survey location(s) to download
    ...				# Additional arguments to RODBC::odbcConnect
){

    if (!requireNamespace("RODBC", quietly = TRUE)) {
        stop("The RODBC package is needed to use this function, but it is not installed. Please install it.", call. = FALSE)
    }  

    start.time <- Sys.time()
    
     params <- c(!missing(loc.prefix), !missing(date.deployed), !missing(date.collected))
    if(sum(params)>0) {
        filter <- " WHERE "
        op <- params
        op[params] <- " AND "
        op[!params] <- ""
    } else {
        filter <- ""
        op <- rep("", 3)
    }
    
    if(!params[1] && all(params[2], params[3])) {
        filter <- " WHERE "
        op[2] <- ""
    }
    
    if(!missing(loc.prefix)){
        if (nchar(loc.prefix[1]) != 6) stop(paste('loc.prefix must be 6 characters, got', loc.prefix))
        loc.prefix <- paste0("(`tblLocation`.`fldLocationNameAbbreviation` = ", paste("'", loc.prefix, "'", sep="", collapse=(" OR `tblLocation`.`fldLocationNameAbbreviation` =  ")), ")")
    } else loc.prefix <- ""
    
    if(!missing(date.deployed)) {
        date.deployed <- paste0("`tblCardRecorder`.`fldDateDeployed` >= '", date.deployed, "'")
    } else date.deployed <- ""
    
    if(!missing(date.collected)) {
        date.collected <- paste0("`tblCardRecorder`.`fldDateCollected` <= '", date.collected, "'")
    } else date.collected <- ""        
    
    query <- paste0("SELECT `tblCardRecorder`.`pkCardRecorderID`, `tblLocation`.`fldLocationNameAbbreviation`, `tblRecorder`.`fldSerialNumber`, `tblCard`.`pkCardID`, `tblCardRecorder`.`fldDateDeployed`, `tblCardRecorder`.`fldDateCollected` FROM `tblCardRecorder` INNER JOIN `tblLocation` ON `tblLocation`.`pkLocationID` = `tblCardRecorder`.`fkLocationID` INNER JOIN `tblRecorder` ON `tblRecorder`.`pkRecorderID` = `tblCardRecorder`.`fkRecorderID` INNER JOIN `tblCard` ON `tblCard`.`pkCardID` = `tblCardRecorder`.`fkCardID`", filter, loc.prefix, op[2], date.deployed, op[3], date.collected, ";")

    # open the database connection
    if(missing(uid) && missing(pwd)) {dbCon <- RODBC::odbcConnect(db.name, ...)
    } else if(missing(uid)) {dbCon <- RODBC::odbcConnect(db.name, pwd, ...)
    } else dbCon <- RODBC::odbcConnect(db.name, uid, pwd, ...)
    
    # Establish a cleanup procedure
    on.exit(close(dbCon))
    
    # Download the query 
    cardrecorder <- RODBC::sqlQuery(dbCon, query)
    
    message(if(class(cardrecorder) == 'data.frame') {paste('Done! Download time:', round(Sys.time()-start.time, 2), 'seconds')
            } else paste("Download unsuccessful; RODBC returned errors: ", paste(cardrecorder, collapse=" ")))
    
    return(cardrecorder)
}    

