# Function to download templates from a database
# Created 2012 Nov 07
# Modified 2015 Sept 6

dbDownloadTemplate <- function(
    db.name='acoustics',                # Connection name in ODBC _and_ on host
    uid,                                # Database User ID, if not in ODBC
    pwd,                                # Database Password, if not in ODBC
    type="BIN",                         # 'BIN', 'binary', or 'bt', also 'COR', 'correlation' or 'ct'
    names,                              # Template names
    species,                            # Template species
    FFTwl,                              # FFT window length for templates
    FFTovlp,                            # FFT overlap
    FFTwn,                              # FFt window name
    ...                                 # Additional arguments to RODBC::odbcConnect
){

    if (!requireNamespace("RODBC", quietly = TRUE)) {
        stop("The RODBC package is needed to use this function, but it is not installed. Please install it.", call. = FALSE)
    }  
    
    start.time <- Sys.time()
    
    # Standarize template type  
    if(tolower(type) %in% c('bin', 'binary', 'bt', 'b')) {std.type <- "'BIN'"
         } else if (tolower(type) %in% c('cor', 'correlation', 'ct', 'c')) {std.type <- "'COR'"
         } else stop(paste('type expects BIN, bt, COR, or ct. Did not recognize:', type))

    params <- c(!missing(names), !missing(species), !missing(FFTwl), !missing(FFTovlp), !missing(FFTwn))
    op <- params
    op[params] <- " AND "
    op[!params] <- ""
    
    if(!missing(names)) names <- paste0("(", paste0("`tblTemplate`.`fldTemplateName` = '", paste0(names, collapse="' OR `tblTemplate`.`fldTemplateName` = '")), "')")
    else names <- ""

    if(!missing(species)) species <- paste0("(", paste0("`tblSpecies`.`fldSpeciesCode` = '", paste0(species, collapse="' OR `tblSpecies`.`fldSpeciesCode` = '")), "')")
    else species <- ""
    
    if(!missing(FFTwl)) FFTwl <- paste0("(", paste0("`tblTemplate`.`fldFFTwl` = '", paste0(FFTwl, collapse="' OR `tblTemplate`.`fldFFTwl` = '")), "')")
    else FFTwl <- ""
    
    if(!missing(FFTovlp)) FFTovlp <- paste0("(", paste0("`tblTemplate`.`fldFFTovlp` = '", paste0(names, collapse="' OR `tblTemplate`.`fldFFTovlp` = '")), "')")
    else FFTovlp <- ""
    
    if(!missing(FFTwn)) FFTwn <- paste0("(", paste0("`tblTemplate`.`fldFFTwn` = '", paste0(FFTwn, collapse="' OR `tblTemplate`.`fldFFTwn` = '")), "')")
    else FFTwn <- ""
    
    
    # Create query
    query <- paste0("SELECT `tblTemplate`.`fldTemplateName`, `tblTemplate`.`fldClipPath`, `tblTemplate`.`fldSampRate`, `tblTemplate`.`fldPtOnT`, `tblTemplate`.`fldPtOnFrq`, `fldPtOffT`, `tblTemplate`.`fldPtOffFrq`, `tblTemplate`.`fldPtsT`, `tblTemplate`.`fldPtsFrq`, `tblTemplate`.`fldPtsAmp`, `tblTemplate`.`fldTStep`, `tblTemplate`.`fldFrqStep`, `tblTemplate`.`fldNTBins`, `tblTemplate`.`fldFirstTBin`, `tblTemplate`.`fldNFrqBins`, `tblTemplate`.`fldDuration`, `tblTemplate`.`fldFrqLim`, `tblTemplate`.`fldFFTwl`, `tblTemplate`.`fldFFTovlp`, `tblTemplate`.`fldFFTwn`, `tblTemplate`.`fldScoreCutoff`, `fldComment` FROM `tblTemplate` INNER JOIN `tblSpecies` ON `tblTemplate`.`fkSpeciesID` = `tblSpecies`.`pkSpeciesID` WHERE `tblTemplate`.`fldActive` = '1' AND `tblTemplate`.`fldTemplateType` = ", std.type, op[1], names, op[2], species, op[3], FFTwl, op[4], FFTovlp, op[5], FFTwn, ";")
    
    # open the database connection
    if(missing(uid) && missing(pwd)) {dbCon <- RODBC::odbcConnect(db.name, ...)
    } else if(missing(uid)) {dbCon <- RODBC::odbcConnect(db.name, pwd, ...)
    } else dbCon <- RODBC::odbcConnect(db.name, uid, pwd, ...)
   
    # Establish a cleanup procedure
    on.exit(close(dbCon))        
    # Download the name and data fields for the templates 
    tblTemplate <- RODBC::sqlQuery(dbCon, query)
    # Assemble pt matrices
    if(std.type == "'BIN'") {
        pt.on <- lapply(1:length(tblTemplate$fldTemplateName), function(x) cbind(t=eval(parse(text=as.character(tblTemplate$fldPtOnT[x]))), frq=eval(parse(text=as.character(tblTemplate$fldPtOnFrq[x])))))
        pt.off=lapply(1:length(tblTemplate$fldTemplateName), function(x) cbind(t=eval(parse(text=as.character(tblTemplate$fldPtOffT[x]))), frq=eval(parse(text=as.character(tblTemplate$fldPtOffFrq[x])))))
    } else if(std.type == "'COR'") {
        pts <- lapply(1:length(tblTemplate$fldTemplateName), function(x) cbind(t=eval(parse(text=as.character(tblTemplate$fldPtsT[x]))), frq=eval(parse(text=as.character(tblTemplate$fldPtsFrq[x]))), amp=-0.01*eval(parse(text=as.character(tblTemplate$fldPtsAmp[x])))))
    }
    
    # reconstruct '__Template' objects
    if(std.type == "'BIN'") {
    	templates <- lapply(1:length(tblTemplate$fldTemplateName), function(x) new('binTemplate', clip.path=as.character(tblTemplate$fldClipPath[x]), samp.rate=tblTemplate$fldSampRate[x], pt.on=pt.on[[x]], pt.off=pt.off[[x]], t.step=tblTemplate$fldTStep[x], frq.step=tblTemplate$fldFrqStep[[x]], n.t.bins=tblTemplate$fldNTBins[x], first.t.bin=tblTemplate$fldFirstTBin[x], n.frq.bins=tblTemplate$fldNFrqBins[x], duration=tblTemplate$fldDuration[x], frq.lim=eval(parse(text=as.character(tblTemplate$fldFrqLim[x]))), wl=tblTemplate$fldFFTwl[x], ovlp=tblTemplate$fldFFTovlp[x], wn=as.character(tblTemplate$fldFFTwn[x]), score.cutoff=tblTemplate$fldScoreCutoff[x], comment=as.character(tblTemplate$fldComment[x])))
    } else if(std.type == "'COR'") {
    	templates <- lapply(1:length(tblTemplate$fldTemplateName), function(x) new('corTemplate', clip.path=as.character(tblTemplate$fldClipPath[x]), samp.rate=tblTemplate$fldSampRate[x], pts=pts[[x]], t.step=tblTemplate$fldTStep[x], frq.step=tblTemplate$fldFrqStep[x], n.t.bins=tblTemplate$fldNTBins[x], first.t.bin=tblTemplate$fldFirstTBin[x], n.frq.bins=tblTemplate$fldNFrqBins[x], duration=tblTemplate$fldDuration[x], frq.lim=eval(parse(text=as.character(tblTemplate$fldFrqLim[x]))), wl=tblTemplate$fldFFTwl[x], ovlp=tblTemplate$fldFFTovlp[x], wn=as.character(tblTemplate$fldFFTwn[x]), score.cutoff=tblTemplate$fldScoreCutoff[x], comment=as.character(tblTemplate$fldComment[x])))
    }        
    # Name each list with the template name    
    names(templates) <- as.character(tblTemplate$fldTemplateName)        
    # Join the templates as '__TemplateList'
    if(std.type == "'BIN'") {
    	templates <- new('binTemplateList', templates=templates)
    } else if(std.type == "'COR'") {
        templates <- new('corTemplateList', templates=templates)
    }
    
    message(if(class(tblTemplate) == 'data.frame') {paste('Done! Download time:', round(Sys.time()-start.time, 2), 'seconds')
            } else paste("Download unsuccessful; RODBC returned errors: ", paste(tblTemplate, collapse=" ")))
    
    return(templates)
}        




