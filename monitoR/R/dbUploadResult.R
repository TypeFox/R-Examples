# Function to upload results to a database
# Created 2012 Oct 23
# Modified 2015 Sept 6

dbUploadResult <- function(
    detection.obj,                   # Pks to upload to database
    which.one,                       # Specify a single element name or number of pks, eg template name, all if blank
    what='detections',               # Upload either 'detections' above the score.cutoff or all 'peaks'
    db.name='acoustics',             # Connection name in ODBC _and_ on host
    uid,                             # Database User ID, if not in ODBC
    pwd,                             # Database Password, if not in ODBC
    analysis.type,                   # Analysis type: COR(relation) or BIN(ary)
    analyst='',                      # From `tblPerson`.`pkPersonID`
    ...                              # Additional arguments to RODBC::odbcConnect
){

    if (!requireNamespace("RODBC", quietly = TRUE)) {
        stop("The RODBC package is needed to use this function, but it is not installed. Please install it.", call. = FALSE)
    }  
    
    start.time <- Sys.time()
    # Pull peaks from 'detectionList' object
    if(missing(which.one)){
        if(tolower(what) %in% c('p', 'peaks', 'pks')) {pks.L <- getPeaks(detection.obj=detection.obj, output='list')
        } else if(tolower(what) %in% c('d', 'detections', 'det')) {pks.L <- getDetections(detection.obj=detection.obj, output='list')
        }
    } else {if(tolower(what) %in% c('p', 'peaks', 'pks')) {pks.L <- getPeaks(detection.obj=detection.obj, which.one=which.one, output='list')
        } else if(tolower(what) %in% c('d', 'detections', 'det')) {pks.L <- getDetections(detection.obj=detection.obj, which.one=which.one, output='list')
        }
    }
    # Assess type
    if(tolower(analysis.type) %in% c('bin', 'bt', 'binary', 'b')) {
    	analysis.type <- 'BIN'
    } else if(tolower(analysis.type) %in% c('cor', 'ct', 'correlation', 'c')) {
    	analysis.type <- 'COR'
    } else 
    	stop('Did not recognize analysis.type, was it BIN or COR?')
    # open the database connection
    if(missing(uid) && missing(pwd)) {dbCon <- RODBC::odbcConnect(db.name, ...)
    } else if(missing(uid)) {dbCon <- RODBC::odbcConnect(db.name, pwd, ...)
    } else dbCon <- RODBC::odbcConnect(db.name, uid, pwd, ...)
    # Establish a cleanup procedure
    on.exit(close(dbCon))
 
    # determine which hits to send to the database
    if (!missing(which.one)) {pks.L <- pks.L[[which.one]]
    } else which.one <- names(pks.L)
    
    # make a new variable to store the score.cutoff for each template
    cutoff <- lapply(detection.obj@templates, function(x) x@score.cutoff)
        
    # make a new variable to record the names of the templates that detected the hits
    templateID <- unlist(lapply(pks.L, function(x) unique(x$template)))
    
    # get list of templates, currently based on fldTemplateName, might change to pkTemplateID??
    template.dat <- RODBC::sqlQuery(dbCon, paste("Select `pkTemplateID`, `fldTemplateName` FROM `tblTemplate` WHERE `fldTemplateName` = '", paste(names(pks.L), sep="", collapse="' OR `fldTemplateName` = '"), "'", sep="")) 
        
    # replace template names in new variable with fkTemplateIDs
    for (i in 1:length(templateID)){
        templateID[i] <- template.dat$pkTemplateID[template.dat$fldTemplateName == templateID[i]]
    }
            
    # get list of surveys, based on fldSurveyName
    survey.name <- detection.obj@survey.name
        
    survey.dat <- RODBC::sqlQuery(dbCon, paste("Select `pkSurveyID`, `fldSurveyName` FROM `tblSurvey` WHERE `fldSurveyName` = '", paste(survey.name, sep="", collapse="' OR `fldSurveyName` = '"), "'", sep="")) 
        
    # replace survey names in new variable with fkSurveyIDs
    for (i in 1:length(survey.name)){
        survey.name[i] <- survey.dat$pkSurveyID[survey.dat$fldSurveyName == survey.name[i]]
    }
        
    # replace pks.L with a list of just the pks or hits
#    pks.L <- lapply(pks.L, '[[', what)
        
    # add a new column containing the templateID to each data frame 
    for (i in 1:length(pks.L)) {pks.L[[i]]$template <- templateID[[i]]}
    
    # add a new column containing the score.cutoff to each data frame
    for (i in 1:length(pks.L)) {pks.L[[i]]$score.cutoff <- cutoff[[i]]} 
        
    # collapse the list into a single data frame
    pks.L <- rbindf(pks.L)
    # convert date.time characters to datetime data type format
    date.time <- unlist(lapply(X=pks.L$date.time, FUN=substr, start=1, stop=19))
    tzone <- unlist(lapply(pks.L$date.time, function(x) as.character(x, format='%Z')))             
    
    # the MySQL query to send the hits to the database
    query<- paste("INSERT INTO `tblResult` (`pkResultID`, `fkSurveyID`, `fkTemplateID`, `fkPersonID`, `fldDateTime`, `fldTimeZone`, `fldTime`, `fldScore`, `fldOnAmp`, `fldOffAmp`, `fldHit`, `fldVerified`, `fldAnalysisType`, `fldLikelihood`, `fldPosterior`, `fldCutoffValue`) VALUES ('", paste(NULL, "', '", survey.name, "', '", pks.L$template, "', '", analyst, "', '", date.time, "', '", tzone, "', '", pks.L$time, "', '", pks.L$score, "', '", if(length(pks.L$on.amp) == 0) {''} else pks.L$on.amp, "', '", if(length(pks.L$off.amp) == 0) {''} else pks.L$off.amp, "', '", if(length(pks.L$hit) == 0) {''} else pks.L$hit*1, "', '", if(length(pks.L$true) == 0) {-1} else pks.L$true*1, "', '", analysis.type, "', '', '', '", pks.L$score.cutoff, "')", sep="", collapse=", ('"), sep="")
    
    # Alert user
    message('\nUploading...')
    # push the query through the open connection to the database     
    status <- RODBC::sqlQuery(dbCon, query)
    
    # report to user
    message(if(is.na(status[1])) {paste('Done! Upload time:', round(Sys.time()-start.time, 2), 'seconds')
            } else if(status[1] == 'character(0)') {paste('Done! Upload time:', round(Sys.time()-start.time, 2), 'seconds')
            } else paste("Upload unsuccessful; RODBC returned errors: ", paste(status, collapse=" ")))
}    


