# Function to download templates from a database
# Modified 2015 Sept 6

dbDownloadResult <- function(
    db.name='acoustics',        # Connection name in ODBC _and_ on host
    uid,                        # Database User ID, if not in ODBC
    pwd,                        # Database Password, if not in ODBC
    survey,                     # File path to survey
    templates,                  # Vector of template names or single templateList object
    type,                       # 'BIN', 'binary', or 'bt', also 'COR', 'correlation' or 'ct'
    FFTwl,                      # optional, for selecting templates
    FFTwn,                      # optional, for selecting templates
    FFTovlp,                    # optional, for selecting templates
    ...                         # Additional arguments to RODBC::sqlQuery
){

    if (!requireNamespace("RODBC", quietly = TRUE)) {
        stop("The RODBC package is needed to use this function, but it is not installed. Please install it.", call. = FALSE)
    }  
    
    start.time <- Sys.time()
      
    # open the database connection
    if(missing(uid) && missing(pwd)) {dbCon <- RODBC::odbcConnect(db.name)
    } else if(missing(pwd)) {dbCon <- RODBC::odbcConnect(db.name, pwd)
    } else dbCon <- RODBC::odbcConnect(db.name, uid, pwd)
    # Establish a cleanup procedure
    on.exit(close(dbCon))
    # Read in survey
    file.ext <- tolower(gsub(".*\\.", "", survey))
    if(file.ext == "wav") {survey.wave <- tuneR::readWave(survey)
    } else if(file.ext == "mp3") survey.wave <- readMP3(survey)
    else stop("survey must be either a wave or mp3 file.")
    # Standarize template type  
    std.type <- if(tolower(type) %in% c('bin', 'binary', 'bt', 'b')) {'BIN'
         } else if (tolower(type) %in% c('cor', 'correlation', 'ct', 'c')) {'COR'
         } else stop(paste('type expects BIN, bt, COR, or ct. Did not recognize:', type))
	# Choose whether to download all templates or just the specified templates
	if(!class(templates) %in% c("binTemplateList", "corTemplateList")) {
            cat('\nDownloading templates\n')
        if(missing(uid) && missing(pwd)) {templates <- dbDownloadTemplate(db.name=db.name, by.cat="names", type=std.type, template.group=templates, FFTwl=FFTwl, FFTwn=FFTwn, FFTovlp=FFTovlp)
        } else if(missing(uid) && !missing(pwd)) {templates <- dbDownloadTemplate(db.name=db.name, pwd=pwd, by.cat="names", type=std.type, template.group=templates, FFTwl=FFTwl, FFTwn=FFTwn, FFTovlp=FFTovlp)
        } else templates <- dbDownloadTemplate(db.name=db.name, uid=uid, pwd=pwd, by.cat="names", type=std.type, template.group=templates, FFTwl=FFTwl, FFTwn=FFTwn, FFTovlp=FFTovlp)
        } 
	# make a new variable to record the names of the templates that detected the hits
	templateID <- templateNames(templates)
	# get list of templates
	template.dat <- RODBC::sqlQuery(dbCon, paste0("Select `pkTemplateID`, `fldTemplateName` FROM `tblTemplate` WHERE `fldTemplateName` = '", paste0(templateID, collapse="' OR `fldTemplateName` = '"), "'")) 
	# replace template names in new variable with fkTemplateIDs
	for (i in 1:length(templateID)){
            templateID[i] <- template.dat$pkTemplateID[template.dat$fldTemplateName == templateID[i]]
	}
	# get list of surveys, based on fldSurveyName
	surveyID <- RODBC::sqlQuery(dbCon, paste0("Select `pkSurveyID`, `fldSurveyName` FROM `tblSurvey` WHERE `fldSurveyName` = '", survey, "'")) 
	# replace survey names in with fkSurveyIDs
	surveyID <- surveyID[surveyID['fldSurveyName'] == survey, 'pkSurveyID']           
    survey.data <- score.L <- peaks.L <- detections.L <- wl <- ovlp <- wn <- list()
    for(i in 1:length(templateID)){
#		score.L[[i]] <- data.frame(
#			date.time=NULL, 
#                    time=NULL, 
#                    score=NULL, 
#                    on.amp=NULL, 
#                    off.amp=NULL
#                    )
		wl[[i]] <- templates@templates[[i]]@wl
		ovlp[[i]] <- templates@templates[[i]]@ovlp
		wn[[i]] <- templates@templates[[i]]@wn   
        if(i == 1) {survey.spec <- spectro(wave=survey.wave, wl=wl[[i]], ovlp=ovlp[[i]], wn=wn[[i]])
        } else if(!all(wl[[i]] == wl[[i-1]], ovlp[[i]] == ovlp[[i-1]], wn[[i]] == wn[[i-1]])) survey.spec <- spectro(wave=survey.wave, wl=wl[[i]], ovlp=ovlp[[i]], wn=wn[[i]])        
        t.bins <- survey.spec$time
        frq.bins <- survey.spec$freq
                        
		survey.data[[i]] <- list(amp=survey.spec$amp, t.bins=t.bins, frq.bins=frq.bins)
		cat('\nDownloading ', templateNames(templates)[i], '\n')
		query <- paste0("SELECT `fldDateTime`, `fldTimeZone`, `fldTime`, `fldOnAmp`, `fldOffAmp`, `fldScore` FROM `tblResult` WHERE `fkSurveyID` = ", surveyID, " AND `fkTemplateID` = ", templateID[i])
#		score.L[[i]] <- RODBC::sqlQuery(dbCon, query, ...)
#		date.time <- paste(score.L[[i]][, 'fldDateTime'], score.L[[i]][, 'fldTimeZone'])
#		score.L[[i]]['fldDateTime'] <- date.time
#		score.L[[i]] <- score.L[[i]][, -2]
#		names(score.L[[i]]) <- c("date.time", "time", "score")
#		peaks.L[[i]] <- score.L[[i]]
        score.L[[i]] <- data.frame("date.time"=NA, "time"=NA, "score"=NA)
		peaks.L[[i]] <- RODBC::sqlQuery(dbCon, query, ...) #
		date.time <- paste(peaks.L[[i]][, 'fldDateTime'], peaks.L[[i]][, 'fldTimeZone']) #
		peaks.L[[i]]['fldDateTime'] <- date.time #
		peaks.L[[i]] <- peaks.L[[i]][, -2] #
		names(peaks.L[[i]]) <- c("date.time", "time", "on.amp", "off.amp", "score") #
		peaks.L[[i]]['detection'] <- FALSE
		peaks.L[[i]][peaks.L[[i]]['score']>=templates@templates[[i]]@score.cutoff, 'detection'] <- TRUE
		detections.L[[i]] <- peaks.L[[i]][peaks.L[[i]]['detection'] == TRUE, c("date.time", "time", "score")]
    }
    names(survey.data) <- templateNames(templates)
    names(score.L) <- names(peaks.L) <- names(detections.L) <- templateNames(templates)
	object <- new('detectionList', survey.name=survey, survey=survey.wave, survey.data=survey.data, templates=templates@templates, scores=score.L, peaks=peaks.L, detections=detections.L)
   cat('\nDone\n')
   return(object)
}        




