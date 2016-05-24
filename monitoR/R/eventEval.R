# Function to compare machine detected events to a verified standard set, and 
# assign an outcome in c("TRUE +", "FALSE +", "TRUE -", "FALSE -") to each event
# Designed to work on only one template at a time; for more than one do a loop
# Written by Jon Katz 14-May 2013

eventEval <- function(
    detections,             # Can be a detection object, data frame, or path to csv file
    what="detections",      # If detections is a detection object specify 'peaks' or 'detections'
    which.one,              # Specify by name or number of which template results to check
    standard,               # File path or data frame of verified detections
    score.cutoff=11,        # Or specify a score cutoff
    tol=1                   # Tolerance of time mismatch; +/- 0.5 * tol seconds from center time
    ){
    
    # cut the tolerance value in half to produce a +/- time buffer
    tol <- 0.5*tol
    # Some minor class checking to extract the detection data
    if(class(detections) == "detectionList") {
        dbDetect <- FALSE
        if(tolower(what) %in% c("detections", "d", "det", "hits")) {getter <- "getDetections"
        } else if(tolower(what) %in% c("peaks", "p", "pks")) getter <- "getPeaks"
        if(missing(which.one) || length(which.one)>1) {stop("Specify 1 'which.one'.")
        } else detections <- do.call(getter, list(detection.obj=detections, which.one=which.one))
        detections <- detections[detections['template'] == which.one, ]
    } else if(class(detections) == "data.frame") {
        if(all(c("fldTime", "fldScore") %in% names(detections))) {dbDetect <- TRUE
        } else if(all(c("time", "score") %in% names(detections))) dbDetect <- FALSE
    } else if(class(detections) == "character") {
        if(tolower(gsub(".*\\.", "", detections)) == 'csv') detections <- read.csv(detections)
        else stop("'detections' must be a detectionList object, a file path to a csv file with a recognizable format, or data frame with a recognizable format.")
        # Check to see if detections were originally from the database
        if("fldTime" %in% names(detections)) {dbDetect <- TRUE
        } else stop('Detections must either be a detectionList object or an equivalent data frame (e.g. from the database).')
    } else if(class(detections) != "data.frame" || all(!c("fldTime", "time", "template") %in% names(detections))) stop("'detections' must be a detectionList object, a file path to a csv file with a recognizeable format, or data frame with a recognizeable format.")
    # Some minor class checking to read in the standard
    if(class(standard) == "character") {detections <- read.csv(standard)
    } else if(class(standard) != "data.frame") stop("'standard' must be a file path or a data frame.")
    # Check to see if standard was originally from the database
    if(all(c("fldTime", "fldTemplateName") %in% names(standard))) {dbStandard <- TRUE
    } else if(all(c("start.time", "name") %in% names(standard))) {dbStandard <- FALSE
    } else stop('Standard must either be from viewSpec() or from the database.')
    # Minimal checking to validate that merging these data frames is even possible...
    # Check to see that all columns exist 
    if(!dbStandard && any(!c("start.time", "end.time", "min.frq", "max.frq", "name") %in% names(standard))) {stop("Unidentified column: ", c("start.time", "end.time", "min.frq", "max.frq", "name")[which(!c("start.time", "end.time", "min.frq", "max.frq", "name") %in% names(standard))], " missing from standard.")
    } else if(dbStandard && !"fldTime" %in% names(standard)) stop("Standard is missing a time field.")
    
    if(nrow(detections)>0 & nrow(standard) == 0) {
        if("fldSurveyName" %in% names(standard) | "fldSurveyName" %in% names(detections)) {
            surv.fld <- "fldSurveyName"
            survey.incl <- TRUE 
            sc.fld <- "fldScore"
        } else if("survey" %in% names(standard) | "survey" %in% names(detections)) {
            surv.fld <- "survey"
            survey.incl <- TRUE
            sc.fld <- "score"
        } else survey.incl <- FALSE
        outcome <- rep("TRUE -", nrow(detections))
        outcome[which(detections[, sc.fld]>=score.cutoff)] <- "FALSE +"
    } else if (any(nrow(detections) == 0, is.null(nrow(detections)))) {
        return(data.frame())
    } else {
        if(dbDetect){ # If detections are from the database
            if(!"fldTime" %in% names(standard)) {
                times <- cbind(standard["start.time"], standard["end.time"])
                times <- rowMeans(times)
                standard['fldTime'] <- times
            }
            if(!"fldSpeciesCode" %in% names(standard)) standard['fldSpeciesCode'] <- standard[, 'name']
            # Build parity between results from the detector and those from the database
            standard[names(detections)[which(!names(detections) %in% names(standard))]] <- NA
            # Reshape the standard data frame
            standard <- standard[names(detections)]
            # Merge the standard with the data
            detections <- rbind(detections, standard)
            # Sort it by time
            detections <- detections[order(detections[, "fldTime"]), ]
            # Determine if the time in each row is the same (+/- tol) as the row above
	        nRows <- nrow(detections)
	        earlytimes <- detections[c(1:(nRows-1)), 'fldTime']
	        latetimes <- detections[c(2:nRows), 'fldTime']
	        manual <- rep(FALSE, nRows)
	        manual[is.na(detections['fldScore'])] <- TRUE
	        manual.above <- c(FALSE, manual[c(1:(nRows-1))])
	        manual.below <- c(manual[c(2:nRows)], FALSE)
            sc.fld <- "fldScore"
            if("fldSurveyName" %in% names(standard) | "fldSurveyName" %in% names(detections)) {
                surv.fld <- "fldSurveyName"
                survey.incl <- TRUE 
            } else survey.incl <- FALSE
        } else if(!dbDetect){ # If detections are straight from the detector
            if(!"time" %in% names(standard)) {
                times <- cbind(standard["start.time"], standard["end.time"])
                times <- rowMeans(times)
                standard['time'] <- times
            }
            if(!"template" %in% names(standard)) standard['template'] <- standard[, 'name']
            # Build parity between results from the detector and those from the database
            standard[names(detections)[which(!names(detections) %in% names(standard))]] <- NA
            # Reshape the standard data frame
            standard <- standard[names(detections)]
            # Merge the standard with the data
            detections <- rbind(detections, standard)
            # Sort it by time
            detections <- detections[order(detections[, "time"]), ]
            # Determine if the time in each row is the same (+/- tol) as the row above
	        nRows <- nrow(detections)
	        earlytimes <- detections[c(1:(nRows-1)), 'time']
	        latetimes <- detections[c(2:nRows), 'time']
	        manual <- rep(FALSE, nRows)
	        manual[is.na(detections['score'])] <- TRUE
	        manual.above <- c(FALSE, manual[c(1:(nRows-1))])
	        manual.below <- c(manual[c(2:nRows)], FALSE)
	        sc.fld <- "score"
            if("survey" %in% names(standard) | "survey" %in% names(detections)) {
                surv.fld <- "survey"
                survey.incl <- TRUE
            } else survey.incl <- FALSE
	    }
	    
	    same1below <- c(earlytimes+tol >= latetimes & earlytimes+tol <= latetimes+tol, FALSE)
	    same1above <- c(FALSE, earlytimes+tol >= latetimes & earlytimes+tol <= latetimes+tol)
	    # Evaluate the time checks and verification checks
	    rows <- as.list(1:nRows)
	    outcome <- unlist(lapply(rows, function(x) {
	        vect <- NULL
	        # The undetected events
		    if(all(manual[x], !same1above[x], !same1below[x])) {vect[x] <- "FALSE -"
		    # Detected events from the standard
		    } else if(manual[x] && any(all(same1above[x], !manual.above[x]), all(same1below[x], !manual.below[x]))) {vect[x] <- NA
		    # Detected but don't exceed score cutoff
		    } else if(!manual[x] && any(all(same1above[x], manual.above[x]), all(same1below[x], manual.below[x])) && detections[x, sc.fld]<score.cutoff) {vect[x] <- "FALSE -"
		    # Detected, exceed score cutoff
		    } else if(!manual[x] && any(all(same1above[x], manual.above[x]), all(same1below[x], manual.below[x])) && detections[x, sc.fld]>=score.cutoff) {vect[x] <- "TRUE +"
		    # Detection exceeding score cutoff not in standard
		    } else if(all(!manual[x], !same1above[x], !same1below[x]) && detections[x, sc.fld]>=score.cutoff) {vect[x] <- "FALSE +"
		    # Detection not exceeding score cutoff & not in standard
		    } else if(all(!manual[x], !same1above[x], !same1below[x]) && detections[x, sc.fld]<score.cutoff) {vect[x] <- "TRUE -"
		    # If above was changed to all(!manual[x], !manual[x-1], !manual[x+1]) would the rest be necessary?
		    } else if(!manual[x] && same1above[x] && !manual.above[x]) {vect[x] <- "TRUE -"
		    } else if(!manual[x] && same1below[x] && !manual.below[x]) {vect[x] <- "TRUE -"
		    } else vect[x] <- NA
		    }))
	}	
		
	# Add outcome to the main data frame
	detections['outcome'] <- outcome
	# Drop the the standard events that were detected 
	detections <- detections[!is.na(detections['outcome']), ]
	# Fill in the survey field for undetected manual events
	if(survey.incl) detections[surv.fld] <- na.omit(unique(detections[, surv.fld]))[1]
	rownames(detections) <- 1:nrow(detections)
    return(detections)
}   

