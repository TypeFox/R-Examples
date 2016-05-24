# Function to upload templates to a database
# Created 2012 Oct 23
# Modified 2015 Sept 6

dbUploadTemplate <- function(
   templates,                           # List of templates to upload to database
   which.one,                           # Can specify one or more template names, or if missing send all of them
   db.name='acoustics',                 # Connection name in ODBC _and_ on host
   uid,                                 # Database User ID, if not in ODBC
   pwd,                                 # Database Password, if not in ODBC
   analyst,                             # From `tblPerson`.`pkPersonID`
   locationID='',                       # from tblsLocation.pkLocationID c(vector)
   date.recorded='',                    # date the song clip was recorded c(vector)
   recording.equip='',                  # recorder (and mic) used to record clip c(vector)
   species.code,                        # vector of sp. codes for templates, e.g. BTNW, ACESAC
   type,                                # types are BIN or COR
   ...                                  # Additional arguments to odbcConnect
){

    if (!requireNamespace("RODBC", quietly = TRUE)) {
        stop("The RODBC package is needed to use this function, but it is not installed. Please install it.", call. = FALSE)
    }  
    
    start.time <- Sys.time()
    # Assess type
    if(tolower(type) %in% c('bin', 'bt', 'binary', 'b')) {
    	type <- 'BIN'
    } else if(tolower(type) %in% c('cor', 'ct', 'correlation', 'c')) {
    	type <- 'COR'
    } else 
    	stop('Did not recognize type, was it BIN or COR?')
    
    # open the database connection
    if(missing(uid) && missing(pwd)) {dbCon <- RODBC::odbcConnect(db.name, ...)
    } else if(missing(pwd)) {dbCon <- RODBC::odbcConnect(db.name, pwd, ...)
    } else dbCon <- RODBC::odbcConnect(db.name, uid, pwd, ...)
    
    # Establish a cleanup procedure
    on.exit(close(dbCon))
    
    # Pull out the templates
    if(missing(which.one)) template.L <- templates@templates
    else template.L <- templates@templates[names(templates@templates) == which.one]
    # Check the species.code vector length
    if(length(species.code)>1 & length(species.code) != length(names(template.L))) stop('You entered ', length(species.code), ' species codes but are uploading ', length(names(template.L)), ' templates, this can\'t be right.')
    # download relevant portions of tblSpecies to lookup fkSpeciesID
    species <- RODBC::sqlQuery(dbCon, paste("SELECT `pkSpeciesID`, `fldSpeciesCode` FROM `tblSpecies` WHERE `fldSpeciesCode` = '", paste(species.code, sep="", collapse="' OR `fldSpeciesCode` = '"), "'", sep=""))
                
    # open a vector to hold fkSpeciesIDs
    speciesID <- NULL
                
    # populate the fkSpeciesIDs
    for (i in 1:length(species.code)){
    speciesID[i] <- species$pkSpeciesID[species$fldSpeciesCode == species.code[i]]
    }
              
    # extract relevant parameters from the templates before upload   
    clips <- lapply(template.L, function(x) x@clip.path)
    srates <- lapply(template.L, function(x) x@samp.rate)
    if(type == 'BIN'){pts.on <- lapply(template.L, function(x) x@pt.on)}
    if(type == 'BIN'){pts.off <- lapply(template.L, function(x) x@pt.off)}
    if(type == 'COR'){pts <- lapply(template.L, function(x) x@pts)}
    t.steps <- lapply(template.L, function(x) x@t.step)
    frq.steps <- lapply(template.L, function(x) x@frq.step)
    n.t.bins <- lapply(template.L, function(x) x@n.t.bins)
    first.t.bin <- lapply(template.L, function(x) x@first.t.bin)
    n.frq.bins <- lapply(template.L, function(x) x@n.frq.bins)
    durations <- lapply(template.L, function(x) x@duration)
    frq.lims <- lapply(template.L, function(x) x@frq.lim)
    wls <- lapply(template.L, function(x) x@wl)
    ovlps <- lapply(template.L, function(x) x@ovlp)
    wns <- lapply(template.L, function(x) x@wn)
    score.cutoffs <- lapply(template.L, function(x) x@score.cutoff)
    comments <- lapply(template.L, function(x) x@comment)
    if(type == 'BIN') {
        pt.on.L <- lapply(template.L, function(x) x@pt.on)
        pt.off.L <- lapply(template.L, function(x) x@pt.off)
        pt.on.t <- lapply(pt.on.L, function(x) x[, 't'])
        pt.on.f <- lapply(pt.on.L, function(x) x[, 'frq'])
        pt.off.t <- lapply(pt.off.L, function(x) x[, 't'])
        pt.off.f <- lapply(pt.off.L, function(x) x[, 'frq'])
    } else if(type == 'COR') {
        pts.L <- lapply(template.L, function(x) x@pts)
        pts.t <- lapply(pts.L, function(x) x[, 't'])
        pts.f <- lapply(pts.L, function(x) x[, 'frq'])
        pts.a <- lapply(pts.L, function(x) x[, 'amp']*-100)
    }

    # the MySQL query to send the template list to the database
    query <- paste0("INSERT INTO `tblTemplate` (`pkTemplateID`, `fkSpeciesID`, `fkPersonID`, `fkLocationID`, `fldTemplateName`, `fldRecordingDate`, `fldRecordingEquipment`, `fldClipPath`, `fldSampRate`, `fldPtOnT`, `fldPtOnFrq`, `fldPtOffT`, `fldPtOffFrq`, `fldPtsT`, `fldPtsFrq`, `fldPtsAmp`, `fldTStep`, `fldFrqStep`, `fldNTBins`, `fldFirstTBin`, `fldNFrqBins`, `fldDuration`, `fldFrqLim`, `fldFFTwl`, `fldFFTovlp`, `fldFFTwn`, `fldScoreCutoff`, `fldTemplateType`, `fldActive`, `fldComment`) VALUES ('", paste0(NULL, "', '", speciesID, "', '", analyst, "', '", locationID, "', '", names(template.L), "', '", date.recorded, "', '", recording.equip, "', '", clips, "', '", srates, "', '", if(type == 'BIN'){pt.on.t} else {''}, "', '", if(type == 'BIN'){pt.on.f} else {''}, "', '", if(type == 'BIN'){pt.off.t} else {''}, "', '", if(type == 'BIN'){pt.off.f} else {''}, "', '", if(type == 'COR'){pts.t} else {''}, "', '", if(type == 'COR'){pts.f} else {''}, "', '", if(type == 'COR'){pts.a} else {''}, "', '", t.steps, "', '", frq.steps, "', '", n.t.bins, "', '", first.t.bin, "', '", n.frq.bins, "', '", durations, "', '", frq.lims, "', '", wls, "', '", ovlps, "', '", wns, "', '", score.cutoffs, "', '", type, "', ", 1,", '", comments, "')", collapse=", ('"))

    # Alert user
    message('Uploading...')      
    # push the query through the open connection to the database     
    status <- RODBC::sqlQuery(dbCon, query)
    # Alert user
    message('Cleaning up...')    
    # Save space in the database and future read vector; lose the white space
    query <- paste("UPDATE `tblTemplate` SET `fldPtOnT` = REPLACE( `fldPtOnT` , ' ' , '' ), `fldPtOnFrq` = REPLACE( `fldPtOnFrq` , ' ' , '' ), `fldPtOffT` = REPLACE( `fldPtOffT` , ' ' , '' ), `fldPtOffFrq` = REPLACE( `fldPtOffFrq` , ' ' , '' ), `fldPtsT` = REPLACE( `fldPtsT` , ' ' , '' ), `fldPtsFrq` = REPLACE( `fldPtsFrq` , ' ' , '' ), `fldPtsAmp` = REPLACE( `fldPtsAmp` , ' ' , '' ) WHERE `fldTemplateName` = '", names(template.L), "'", sep="")
    lapply(query, function(x) RODBC::sqlQuery(dbCon, x))
    # And drop the newline too
    query <- paste("UPDATE `tblTemplate` SET `fldPtOnT` = REPLACE( `fldPtOnT` , '\n' , '' ), `fldPtOnFrq` = REPLACE( `fldPtOnFrq` , '\n' , '' ), `fldPtOffT` = REPLACE( `fldPtOffT` , '\n' , '' ), `fldPtOffFrq` = REPLACE( `fldPtOffFrq` , '\n' , '' ), `fldPtsT` = REPLACE( `fldPtsT` , '\n' , '' ), `fldPtsFrq` = REPLACE( `fldPtsFrq` , '\n' , '' ), `fldPtsAmp` = REPLACE( `fldPtsAmp` , '\n' , '' ) WHERE `fldTemplateName` = '", names(template.L), "'", sep="")
    lapply(query, function(x) RODBC::sqlQuery(dbCon, x))
    # Might as well drop carriage returns too
    query <- paste("UPDATE `tblTemplate` SET `fldPtOnT` = REPLACE( `fldPtOnT` , '\r' , '' ), `fldPtOnFrq` = REPLACE( `fldPtOnFrq` , '\r' , '' ), `fldPtOffT` = REPLACE( `fldPtOffT` , '\r' , '' ), `fldPtOffFrq` = REPLACE( `fldPtOffFrq` , '\r' , '' ), `fldPtsT` = REPLACE( `fldPtsT` , '\r' , '' ), `fldPtsFrq` = REPLACE( `fldPtsFrq` , '\r' , '' ), `fldPtsAmp` = REPLACE( `fldPtsAmp` , '\r' , '' ) WHERE `fldTemplateName` = '", names(template.L), "'", sep="")
    lapply(query, function(x) RODBC::sqlQuery(dbCon, x))
    
    # report to user
    message(if(is.na(status[1])) {paste('Done! Upload time:', round(Sys.time()-start.time, 2), 'seconds')
            } else if(status[1] == 'character(0)') {paste('Done! Upload time:', round(Sys.time()-start.time, 2), 'seconds')
            } else paste("Upload unsuccessful; RODBC returned errors: ", paste(status, collapse=" ")))
}




