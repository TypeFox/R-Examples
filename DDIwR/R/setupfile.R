setupfile <- function(lbls = "", type="all", csv = "", miss, trymiss = FALSE, uniqueid = "",
                      SD = "", delimiter = ",", OS = "windows", outfile="", ...) {
    
    # change the intrnlbls argument into a very unique name, just in case the
    # list object in the external R file(s) are named "intrnlbls" as well
    
    intrnlbls <- lbls
    intrnlbls_objname <- deparse(substitute(lbls))
    rm(lbls)
    
    other.args <- list(...)
    
    pathIsFolder <- FALSE
    if ("pathIsFolder" %in% names(other.args)) {
        pathIsFolder <- other.args$pathIsFolder
    }
    
    saveFile <- FALSE
    if ("saveFile" %in% names(other.args)) {
        saveFile <- other.args$saveFile
    }
    
    
    if (OS == "") {
        OS <- Sys.info()[['sysname']]
    }
    
    if (all(is.character(intrnlbls))) { # all() just in case someone provides a vector by mistake
        
        if (length(intrnlbls) > 1) {
            cat("\n")
            stop("The lbls argument should contain a single path to the list object.\n\n", call. = FALSE)
        }
        
        xmlfiles <- FALSE
        
        labelist <- treatPath(intrnlbls, type = "R")
        
        if (length(labelist) == 1) {
            labelist <- treatPath(intrnlbls, type = "XML")
            if (length(labelist) == 1) {
                cat("\n")
                stop(gsub("XML", "R or .XML", labelist), call. = FALSE)
            }
            else {
                xmlfiles <- TRUE
            }
        }
        
        if (!file.exists("Setup files")) {
            dir.create("Setup files")
        }
        
        csvdatadir <- FALSE # by default
        
        
        # now trying to assess what the csv argument is
        # it can be an object containing csv data, or
        # it can be a string containing a path to the data
        
        
        
        if (all(is.character(csv))) {
            if (csv != "") {
                if (length(csv) > 1) {
                    cat("\n")
                    stop("The csv argument should contain a single path to the list object.\n\n", call. = FALSE)
                }
                
                csvlist <- treatPath(csv, type = "csv")
                if (length(csvlist) > 1) {
                    datadir <- csvlist$completePath
                    csvdatadir <- TRUE
                }
                else {
                    cat("\nNOTE:", csvlist)
                          # since "csvlist" is now an error message from treatPath()
                }
            }
            else {
                # it's important to differentiate between "data" and "Data", for OSs that are case sensitive
                csvdatadir <- file.exists(datadir <- file.path(labelist$completePath, "data"))
                datathere <- csvdatadir
                csvdatadir <- file.exists(datadir <- file.path(labelist$completePath, "Data"))
                
                if (csvdatadir) {
                    csvlist <- treatPath(datadir, type = "csv")
                    if (length(csvlist) == 1) {
                        csvdatadir <- FALSE
                        cat(paste("\nNOTE: There is a ", ifelse(datathere, "data", "Data"), " directory within ",
                            labelist$completePath, ". "), csvlist)
                    }
                }
            }
        }
        
        
        if (csvdatadir) {
            csvfiles <- csvlist$files
            csvnames <- csvlist$filenames
            csvext <- csvlist$fileext
            
            cat ("Processing (including data directory):\n")
        }
        else {
            
            if (is.data.frame(csv)) {
                if (length(labelist$files) > 1) {
                    cat("\n")
                    stop("There are multiple files containing labels and only one csv file provided.\n\n", call. = FALSE)
                }
            }
            
            cat("Processing (no data directory):\n")
        }
        
        
        for (i in seq(length(labelist$files))) {
            
            if (xmlfiles) {
                intrnlblsObject <- getMetadata(file.path(labelist$completePath, labelist$files[i]), fromsetupfile = TRUE, saveFile = saveFile)
            }
            else {
                aa <- ls()
                
                tryCatch(eval(parse(file.path(labelist$completePath, labelist$files[i]))), error = function(x) {
                    stop(paste("\nThere is an error associated with the file \"", labelist$files[i], "\", see below:\n       ", gsub("Error in ", "", as.character(x)), sep=""), call. = FALSE)
                })
                
                bb <- ls()
                bb <- bb[-which(bb == "aa")]
            }
            
            if (csvdatadir) {
                
                if (labelist$filenames[i] %in% csvnames) {
                    cat(labelist$filenames[i], "\n")
                    position <- match(labelist$filenames[i], csvnames)
                    
                    
                    for (j in seq(length(position))) {
                        
                        if (csvext[position[j]] %in% c("CSV", "CSV.GZ")) {
                                                                                                   # delimiter is already set from the function's formal argument
                            csvreadfile <- read.csv(file.path(datadir, csvfiles[position[j]]), sep = delimiter, header = TRUE, as.is = TRUE)
                            
                            if (ncol(csvreadfile) == 1) {
                                delimiter <- getDelimiter(file.path(datadir, csvfiles[position[j]]))
                                
                                if (delimiter == "unknown") {
                                    stop(paste("Unknown column separator for the file", csvfiles[position[j]],
                                               "\nShould be either \",\" or \";\" or tab separated.\n\n"), call. = FALSE)
                                }
                                
                                csvreadfile <- read.csv(file.path(datadir, csvfiles[position[j]]), sep = delimiter, header = TRUE, as.is = TRUE)
                            }
                            
                            
                            if (!xmlfiles) {
                                intrnlblsObject <- get(setdiff(bb, aa))
                            }
                            tryCatch(Recall(intrnlblsObject, type = type, miss = miss, csv = csvreadfile, trymiss = trymiss, uniqueid = uniqueid, SD = SD,
                                            delimiter = delimiter, OS = OS, outfile = labelist$filenames[i], pathIsFolder = pathIsFolder, ... = ...),
                                error = function(x) {
                                    cat(paste("     There is an error associated with the file \"", labelist$filenames[i], "\", see below:\n     ", sep=""))
                                    cat(as.character(x))
                                })
                        }
                    }
                }
                else {
                    cat(labelist$filenames[i], "(no .csv file)", "\n")
                    if (!xmlfiles) {
                        intrnlblsObject <- get(setdiff(bb, aa))
                    }
                    tryCatch(Recall(intrnlblsObject, type = type, miss = miss, trymiss = trymiss, uniqueid = uniqueid, SD = SD, 
                                    delimiter = delimiter, OS = OS, outfile = labelist$filenames[i], pathIsFolder = pathIsFolder, ... = ...),
                        error = function(x) {
                            cat(paste("     There is an error associated with the file \"", labelist$filenames[i], "\", see below:\n     ", sep=""))
                            cat(as.character(x))
                        })
                }
            }
            else {
                cat(labelist$filenames[i], "\n")
                
                if (is.data.frame(csv)) {
                    if (length(labelist$filenames) == 1) {
                        if (!xmlfiles) {
                            intrnlblsObject <- get(setdiff(bb, aa))
                        }
                        tryCatch(Recall(intrnlblsObject, type = type, miss = miss, csv = csv, trymiss = trymiss, uniqueid = uniqueid, SD = SD,
                                        delimiter = delimiter, OS = OS, outfile = labelist$filenames[i], pathIsFolder = pathIsFolder, ... = ...),
                        error = function(x) {
                            cat(paste("     There is an error associated with the file \"", labelist$filenames[i], "\", see below:\n     ", sep=""))
                            cat(as.character(x))
                        })
                    }
                }
                else {
                    # there is really no csv data
                        if (!xmlfiles) {
                            intrnlblsObject <- get(setdiff(bb, aa))
                        }
                        tryCatch(Recall(intrnlblsObject, type = type, miss = miss, trymiss = trymiss, uniqueid = uniqueid, SD = SD, 
                                        delimiter = delimiter, OS = OS, outfile = labelist$filenames[i], pathIsFolder = pathIsFolder, ... = ...),
                        error = function(x) {
                            cat(paste("     There is an error associated with the file \"", labelist$filenames[i], "\", see below:\n     ", sep=""))
                            cat(as.character(x))
                        })
                }
            }
            
            if (!xmlfiles) {
                rm(list = c(eval(setdiff(bb, aa)), "bb", "aa"))
            }
        }
        
        cat("\nSetup files created in:\n", file.path(getwd(), "Setup files"), "\n\n", sep="")
        
        return(invisible())
    }
    
    
    
    
    csvlist <- NULL # initialization
    if (all(is.character(csv))) {
        if (all(csv != "")) {
            if (length(csv) > 1) {
                cat("\n")
                stop("The csv argument should contain a single path to the list object.\n\n", call. = FALSE)
            }
            
            csvlist <- treatPath(csv, type = "CSV")
            if (length(csvlist) > 1) {
                # no error
                if (length(csvlist$files) > 1) {
                    cat("\n")
                    stop("There is only one object containing labels and multiple csv files.\n\n", call. = FALSE)
                }
            }
            else {
                # There is a single string returned by treatPath(), with an error message
                cat("\nNOTE:", csvlist)
                csv <- "" # back to the default value
            }
        }
    }
    
    
    
    if (is.null(names(intrnlbls)) | !all(names(intrnlbls) %in% c("varlab", "vallab"))) {
        cat("\n")
        stop("The object does not contain labels for variables and/or values.\n\n", call. = FALSE)
    }
    
    if (!(type %in% c("SPSS", "Stata", "SAS", "R", "all"))) {
        cat("\n")
        stop("The argument <type> can only be: \"SPSS\", \"Stata\", \"SAS\", \"R\", or \"all\".\n\n", call. = FALSE)
    }
    
    enter <- getEnter(OS=OS)
    
    varnames <- names(intrnlbls$varlab)
    maxchars <- max(nchar(varnames))
    varcheck <- rep(0, length(varnames))
    formats <- FALSE
    
    csv_is_df <- is.data.frame(csv)
    
    csv_is_path <- FALSE
    if (length(csv) == 1) { # csv is a character vector of length 1, i.e. a path 
        if (is.character(csv)) {
            if (csv != "") {
                csv_is_path <- TRUE
            }
        }
    }
    
    if (csv_is_df | csv_is_path) {
        
        if (!is.null(csvlist)) {
                                                                                           # delimiter is already set from the function's formal argument
            csvreadfile <- read.csv(file.path(csvlist$completePath, csvlist$files[1]), sep = delimiter, header = TRUE, as.is=TRUE)
            
            if (ncol(csvreadfile) == 1) {
            
                delimiter <- getDelimiter(file.path(csvlist$completePath, csvlist$files[1]))
                
                if (delimiter == "unknown") {
                    stop(paste("Unknown column separator for the file", csvlist$files[1],
                               "\nShould be either \",\" or \";\" or tab separated.\n\n"), call. = FALSE)
                }
                
                csvreadfile <- read.csv(file.path(csvlist$completePath, csvlist$files[1]), sep = delimiter, header = TRUE, as.is=TRUE)
            }
            
            # cat("\n")
            # cat("Found \"", csvlist$files[1], "\" in the directory \"", csvlist$completePath, "\". Using that as the .csv file.\n\n", sep="")
            
            csv <- csvreadfile
        }
        
        csvnames <- names(csv)
        csvformats <- sasformats <- rep("", length(csvnames))
        if (!is.data.frame(csv)) {
            cat("\n")
            stop("The csv file should be a data frame.\n\n", call. = FALSE)
        }
        
        gofurther <- TRUE
        
        plusnames <- setdiff(toupper(csvnames), toupper(names(intrnlbls$varlab)))
        if (length(plusnames) > 0) {
            if (length(plusnames) == length(csvnames)) {
                cat("    None of the variables in the .csv file have metadata information.\n",
                    "    (perhaps the .csv file doesn't have the variable names in the first row?)\n", sep="")
                gofurther <- FALSE
            }
            else {
                cat("    There is no metadata information for the following variables in the .csv file:\n")
                plusnames <- strwrap(paste(plusnames, collapse=", "), 75)
                for (pnms in plusnames) {
                    cat("       ", pnms, "\n")
                }
                cat("\n")
            }
        }
        
        
        plusnames <- setdiff(toupper(names(intrnlbls$varlab)), toupper(csvnames))
        if (length(plusnames) > 0) {
            cat("    There is metadata information for the following variables, but *not* in the .csv file:\n")
            plusnames <- strwrap(paste(plusnames, collapse=", "), 75)
            for (pnms in plusnames) {
                cat("       ", pnms, "\n")
            }
            
            if (gofurther) {
                cat("       ", ifelse(length(plusnames) == 1, "This variable", "These variables"), "will be omitted.\n")
            }
            else {
                cat("\n")
            }
        }
        
        nrowscsv <- nrow(csv)
        
        
        
        if (gofurther) {
            
            intrnlbls$varlab <- intrnlbls$varlab[which(toupper(names(intrnlbls$varlab)) %in% toupper(csvnames))]
            intrnlbls$vallab <- intrnlbls$vallab[which(toupper(names(intrnlbls$vallab)) %in% toupper(csvnames))]
            
            varnames <- names(intrnlbls$varlab)
            maxchars <- max(nchar(varnames))
            varcheck <- rep(0, length(varnames))
            
            printNOTE <- FALSE
            
            for (i in seq(length(csvnames))) {
                vartype <- "numeric"
                decimals <- FALSE
                
                tempvar <- csv[, csvnames[i]]
                
                if (is.character(tempvar)) {
                    vartype <- "string"
                    if (any(tempvar == ".")) { # Stata type empty cells
                        tempvar[tempvar == "."] <- NA
                        printNOTE <- TRUE
                    }
                }
                
                nofchars <- nchar(as.character(tempvar))
                nofchars[is.na(tempvar)] <- 0
                maxvarchar <- max(nofchars)
                
                if (toupper(csvnames[i]) %in% names(intrnlbls$vallab)) {
                    if (is.character(intrnlbls$vallab[[toupper(csvnames[i])]])) {
                        vartype <- "string"
                        # gofurther <- FALSE
                        maxvarchar <- max(maxvarchar, nchar(intrnlbls$vallab[[toupper(csvnames[i])]]))
                    }
                }
                
                if (all(is.na(tempvar))) { # completely empty variable
                    vartype <- "missing"
                }
                
                if (vartype == "numeric") {
                    if (max(tempvar, na.rm = TRUE) - floor(max(tempvar, na.rm = TRUE)) > 0) { # has decimals
                        decimals <- TRUE
                    }
                    else {
                        if (toupper(csvnames[i]) %in% names(intrnlbls$vallab)) {
                            if (is.numeric(intrnlbls$vallab[[toupper(csvnames[i])]])) {
                                maxvarchar <- max(maxvarchar, nchar(length(intrnlbls$vallab[[toupper(csvnames[i])]])))
                            }
                        }
                    }
                    
                    if (decimals) {
                        csvformats[i] <- paste("F", maxvarchar, ".2", sep="")
                    }
                    else {
                        csvformats[i] <- paste("F", maxvarchar, ".0", sep="")
                    }
                }
                else if (vartype == "string") {
                    sasformats[i] <- "$"
                    csvformats[i] <- paste("A", maxvarchar, sep="")
                }
                else { # all the values are missing
                    csvformats[i] <- paste("F1.0", sep="")
                }
                
                varcheck[i] <- 1
            }
            
            if (printNOTE) {
                cat("    NOTE: some variable(s) in this file have a \".\" sign to represent a missing.\n")
                cat("    The .csv file might not import properly in some software.\n\n")
            }
            
            formats <- all(csvformats != "")
        }
        
        ## TO check if any of the existing metadata variables is not found in the CSV data file
    }
    
    stringvars <- lapply(intrnlbls$vallab, function(x) {
        all(is.character(x[[1]]))
    })
    
    
    if (missing(outfile)) {
        if (grepl("\"", intrnlbls_objname)) {
            outfile <- readline("Name for the setup file:\n")
        }
        else {
            outfile <- intrnlbls_objname
        }
    }
    
    uniqueList <- lapply(unique(intrnlbls$vallab), function(uniques) {
        vars <- sapply(names(intrnlbls$vallab),
                     function(x) {
                         ifelse(length(intrnlbls$vallab[[x]]) == length(uniques),
                                all(names(intrnlbls$vallab[[x]]) == names(uniques)), FALSE)
                     })
        return(names(vars[vars]))
    })
    
    
    if (missing(miss)) {
        if (trymiss) {
            miss <- c("DK/NA", "DK/NO", "DK", "NA", "N/A", "N.A.", "Not answered",
                      "Don't know", "(Don't know)", "No answer", "No opinion",
                      "Not applicable", "Not relevant", "Refused", "(Refused)",
                      "Refused / no answer", "(Refused / no answer)",
                      "Can't say", "Don't know / Can't say")
        }
    }
    
    
    if (type == "SPSS" | type == "all") {
        intrnlbls2 <- intrnlbls
        printMISSING <- FALSE
        
        if (!file.exists("Setup files")) {
            dir.create("Setup files")
        }
        
        if (!file.exists(file.path("Setup files", "SPSS"))) {
            dir.create(file.path("Setup files", "SPSS"))
        }
        
        
        currentdir <- getwd()
        setwd(file.path("Setup files", "SPSS"))
        sink(ifelse(length(grep("\\.sps", outfile)) > 0, outfile, paste(outfile, ".sps", sep="")))
        
        cat("* ------------------------------------------------------------------------------", enter, enter,
            "* --- CONFIGURATION SECTION - START ---", enter, enter, enter, sep="")

        if (formats) {
            cat("* The following command should contain the complete path and", enter,
                "* name of the .csv file to be read (e.g. \"C:/CIS 2008/Data/ALL.csv\")", enter,
                "* Change CSV_DATA_PATH to your filename, below:", enter, enter,
                "FILE HANDLE csvpath /NAME=\"CSV_DATA_PATH\" .", enter, enter, enter, sep="")
        }
        
        cat("* The following command should contain the complete path and", enter,
            "* name of the .sav file to be saved (e.g. \"C:/CIS 2008/Data/ALL.sav\")", enter,
            "* Change SAV_DATA_PATH to your filename, below:", enter, enter,
            "FILE HANDLE savfile /NAME=\"SAV_DATA_PATH\" .", enter, enter, enter,
            "* --- CONFIGURATION SECTION -  END  ---", enter, enter,
            "* ------------------------------------------------------------------------------", enter, enter, enter, enter,
            "* There should be nothing to change below this line", enter,                                                            
            "* ------------------------------------------------------------------------------", enter, enter, enter, enter, sep="")
                  
        if (formats) {
            cat("* -------------- Start Definition Macro --------------", enter, enter,
                "SET LOCALE = 'English' .", enter,
                "SHOW LOCALE .", enter, enter, # SET DECIMAL = COMMA . * (might be another idea)
                "* --------------     Read Raw Data      --------------", enter, enter,
                "GET DATA", enter,
                " /TYPE=TXT", enter,
                " /FILE=csvpath", enter,
                " /DELCASE=LINE", enter,
                " /DELIMITERS=\"", ifelse(delimiter == "\t", "\\t", delimiter), "\"", enter,
                " /ARRANGEMENT=DELIMITED", enter,
                " /FIRSTCASE=2", enter,
                " /IMPORTCASE=ALL", enter,
                " /VARIABLES=", enter, sep="")
            
            maxcharcsv <- max(nchar(csvnames))
            for (i in seq(length(csvnames))) {
                cat(toupper(csvnames[i]), paste(rep(" ", maxcharcsv - nchar(csvnames[i]) + 1), collapse=""), csvformats[i], sep="")
                if (i == length(csvnames)) {
                    cat(" .")
                }
                cat(enter)
            }
            cat("CACHE .", enter, "EXECUTE .", enter, enter,
                "* ------------------------------------------------------------------------------", enter, enter, enter, sep="")
        }
        
        
        if (any(unlist(stringvars))) {
            cat("* --- Recode string variables which have labels, to numeric variables ---", enter, enter, enter)
            for (i in names(stringvars)) {
                if (stringvars[[i]]) {
                    
                    precommand <- paste("RECODE ", toupper(i), " ", sep="")
                    postcommand <- "(\""
                    
                    command <- paste(precommand, postcommand, sep = "")
                    vals <- seq(length(intrnlbls2$vallab[[i]]))
                    for (j in vals) {
                        postcommand <- paste(postcommand, intrnlbls2$vallab[[i]][j], "\"", " = ", j, sep = "")
                        command <- paste(command, intrnlbls2$vallab[[i]][j], "\"", " = ", j, sep = "")
                        if (j == length(intrnlbls2$vallab[[i]])) {
                            command <- paste(command, ") INTO TEMPVRBL .", enter, "EXECUTE .", enter, enter, sep="")
                        }
                        else {
                            if (nchar(postcommand) > 70) {
                                postcommand <- paste(paste(rep(" ", nchar(precommand)), collapse=""), "(\"", sep="")
                                command <- paste(command, ")", enter, paste(rep(" ", nchar(precommand)), collapse=""), "(\"", sep="")
                            }
                            else {
                                command <- paste(command, ") (\"", sep="")
                            }
                        }
                    }
                    
                    cat(command)
                    
                    cat("ALTER TYPE ", toupper(i), "(F", max(nchar(vals)),".0) .", enter, enter,
                        "COMPUTE ", toupper(i), " = TEMPVRBL .", enter, "EXECUTE .", enter, enter, sep="")
                    cat("DELETE VARIABLES TEMPVRBL .", enter, "EXECUTE .", enter, enter, sep="")
                    
                    names(vals) <- names(sort(intrnlbls2$vallab[[i]]))
                    intrnlbls2$vallab[[i]] <- vals
                }
            }
            cat(enter)
        }
        
        cat("* --- Add variable labels --- ", enter, enter,
            "VARIABLE LABELS", enter, sep="")
        
        for (i in seq(length(varnames))) {
            cat(toupper(varnames[i]), paste(rep(" ", maxchars - nchar(varnames[i])), collapse=""), " \"", intrnlbls2$varlab[[i]][1], "\"", sep="")
            if (i == length(varnames)) {
                cat(" .", enter, "EXECUTE .", sep="")
            }
            cat(enter)
        }
        cat(enter, enter, sep="")
        
        cat("* --- Add value labels --- ", enter, enter,
            "VALUE LABELS", enter, sep="")
    
        for (i in seq(length(uniqueList))) {
            n <- uniqueList[[i]][1]
            
            cat(splitrows(uniqueList[[i]], enter, 80), enter, sep="")
            
            #if (all(is.character(intrnlbls2$vallab[[n]]))) {
            #    cat(paste(paste("\"", intrnlbls2$vallab[[n]], "\" \"", names(intrnlbls2$vallab[[n]]), "\"", sep=""), collapse="\n"))
            #}
            #else {
                cat(paste(paste(intrnlbls2$vallab[[n]], " \"", names(intrnlbls2$vallab[[n]]), "\"", sep=""), collapse=enter))
            #}
            if (i == length(uniqueList)) {
                cat(" .", enter, "EXECUTE .", enter, sep="")
            }
            else {
                cat(enter, "/", enter, sep="")
            }
        }
        cat(enter, enter, sep="")
        
        
        if (!missing(miss)) {
        
            if (is.numeric(miss)) {
                missvars <- lapply(intrnlbls2$vallab, function(x) !is.na(match(miss, x)))
                withmiss <- as.vector(unlist(lapply(missvars, any)))
                missvals <- unique(lapply(missvars[withmiss], function(x) miss[x]))
            }
            else {
                missvars <- lapply(intrnlbls2$vallab, function(x) !is.na(match(names(x), miss)))
                withmiss <- as.vector(unlist(lapply(missvars, any)))
                missvals <- unique(lapply(which(withmiss), function(x) as.vector(intrnlbls2$vallab[[x]][missvars[[x]]])))
            }
            
            msngs <- unique(unlist(missvals))
            
            withmiss2 <- which(withmiss)
            
            if(length(missvals) > 0) {
                uniqueMissList <- list()
                for (i in seq(length(missvals))) {
                    vars <- NULL
                    for (j in withmiss2) {
                        y <- intrnlbls2$vallab[[j]][missvars[[j]]]
                        if (all(y %in% missvals[[i]])) {
                            vars <- c(vars, names(intrnlbls2$vallab)[j])
                        }
                    }
                    uniqueMissList[[i]] <- vars
                }
                
                cat("* --- Add missing values --- ", enter, enter,
                    "MISSING VALUES", enter, sep="")
                    
                for (i in seq(length(uniqueMissList))) {
                    if (length(missvals[[i]]) < 4) {
                        cat(splitrows(uniqueMissList[[i]], enter, 80))
                        cat(" (", paste(missvals[[i]], collapse=", ") , ")", sep="")
                    }
                    else {
                        absrange <- abs(range(missvals[[i]]))
                        
                        if (all(missvals[[i]] < 0)) {
                            cat(splitrows(uniqueMissList[[i]], enter, 80))
                            cat(" (LOWEST THRU ", max(missvals[[i]]) , ")", sep="")
                        }
                        else {
                            # check if the missing values range doesn't contain any other (non-missing) values
                            checklist <- list()
                            for (mv in uniqueMissList[[i]]) {
                                allvalues <- intrnlbls2$vallab[[mv]]
                                nonmiss <- allvalues[!allvalues %in% missvals[[i]]]
                                checklist[[mv]] <- any(nonmiss %in% seq(min(missvals[[i]]), max(missvals[[i]])))
                            }
                            
                            checklist <- unlist(checklist)
                            
                            # print(checklist)
                            
                            if (any(checklist)) {
                                # at least one variable has a non-missing value within the range of the missing values
                                printMISSING <- TRUE
                                
                                # now trying to see if at least some of the variables can be "rescued"
                                if (any(!checklist)) {
                                    ###
                                    # poate merge...? DE TESTAT
                                    cat(splitrows(names(checklist)[!checklist], enter, 80))
                                    # cat(paste(names(checklist)[!checklist], collapse=", "))
                                    ###
                                    cat(paste(" (", min(missvals[[i]]), " TO ", max(missvals[[i]]) , ")", enter, sep=""))
                                    checklist <- checklist[checklist]
                                }
                                
                                ###
                                # poate merge...? DE TESTAT
                                cat(splitrows(names(checklist), enter, 80))
                                # cat(paste(names(checklist), collapse=", "))
                                ###
                                cat(paste(" (", paste(missvals[[i]][1:3], collapse=", ") , ")", sep=""))
                                cat(ifelse(i == length(uniqueMissList), " .", ""))
                                cat("  * more than three distinct missing values found")
                            }
                            else {
                                cat(splitrows(uniqueMissList[[i]], enter, 80))
                                cat(" (", min(missvals[[i]]), " TO ", max(missvals[[i]]) , ")", sep="") 
                            }
                        }
                    }
                    cat(ifelse(i == length(uniqueMissList), " .", ""))
                    cat(enter)
                }
                cat(enter, enter, sep="")
            }
        }
        
        cat(paste("* --- Save the .sav file --- ", enter, enter,
                  "SAVE OUTFILE=savfile", enter, "/KEEP", enter, sep=""))
        
        if (formats) {
            for (n in csvnames) {
                cat(toupper(n), enter, sep="")
            }
        }
        else {
            for (n in names(intrnlbls2$varlab)) {
                cat(toupper(n), enter, sep="")
            }
        }
        
        cat("  /COMPRESSED .", enter, "EXECUTE .", enter, sep="")
        
        # finish writing and close the .sps file
        sink()
        
        if (printMISSING) {
            # this would be printed on the screen
            cat("    For some variables, more than 3 distinct missing values were found.\n")
            cat("    Only the first three were used.\n\n")
        }
        
        setwd(currentdir)
        
    }
    
    
    
    if (type == "Stata" | type == "all") {
        intrnlbls2 <- intrnlbls
        uniqueList2 <- uniqueList
        
        if (!file.exists("Setup files")) {
            dir.create("Setup files")
        }
        
        if (!file.exists(file.path("Setup files", "Stata"))) {
            dir.create(file.path("Setup files", "Stata"))
        }
        
        currentdir <- getwd()
        setwd(file.path("Setup files", "Stata"))
        sink(ifelse(length(grep("\\.do", outfile)) > 0, outfile, paste(outfile, ".do", sep="")))
        
        cat("/* Initialization commands */", enter,
            "clear", enter,
            "capture log close", enter,
            "set more off", enter,
            "version 12.0", enter,
            "set linesize 250", enter,
            "set varabbrev off", enter,
            "set memory 1G // not necessary in Stata 12",
            ifelse(SD == ";", "\n#delimit ;", ""), enter, enter, enter,
            "* ----------------------------------------------------------------------------",
            ifelse(SD == ";", " ;", ""), enter, enter,
            "* --- CONFIGURATION SECTION - START ---",
            ifelse(SD == ";", "                                        ;", ""), enter, enter,
            "* The following command should contain the complete path and",
            ifelse(SD == ";", "                   ;", ""), enter,
            "* name of the Stata log file.",
            ifelse(SD == ";", "                                                  ;", ""), enter,
            "* Change LOG_FILENAME to your filename, below:", enter,
            ifelse(SD == ";", "                                 ;", ""), enter,
            "local log_file \"LOG_FILENAME\"",
            ifelse(SD == ";", " ;", ""), enter, enter, enter,
            "* The following command should contain the complete path and",
            ifelse(SD == ";", "                   ;", ""), enter,
            "* name of the CSV file, usual file extension \".csv\"",
            ifelse(SD == ";", "                            ;", ""), enter,
            "* Change CSV_DATA_PATH to your filename, below:", enter,
            ifelse(SD == ";", "                                 ;", ""), enter,
            "local csvpath \"CSV_DATA_PATH\"",
            ifelse(SD == ";", " ;", ""), enter, enter, enter,
            "* The following command should contain the complete path and",
            ifelse(SD == ";", "                   ;", ""), enter,
            "* name of the STATA file, usual file extension \".dta\"",
            ifelse(SD == ";", "                          ;", ""), enter,
            "* Change STATA_DATA_PATH to your filename, below:", enter,
            ifelse(SD == ";", "                               ;", ""), enter,
            "local statapath \"STATA_DATA_PATH\"",
            ifelse(SD == ";", " ;", ""), enter, enter, enter,
            "* --- CONFIGURATION SECTION - END ---",
            ifelse(SD == ";", "                                          ;", ""), enter, enter,
            "* ----------------------------------------------------------------------------",
            ifelse(SD == ";", " ;", ""), enter, enter, enter, enter,
            "* There should be nothing to change below this line",
            ifelse(SD == ";", "                   ;", ""), enter,
            "* ----------------------------------------------------------------------------",
            ifelse(SD == ";", " ;", ""), enter, enter, enter, enter,         
            "log using \"`log_file'\", replace text",
            ifelse(SD == ";", " ;", ""), enter, enter,
            "insheet using \"`csvpath'\", comma names case",
            ifelse(SD == ";", " ;", ""), enter, enter,
            "* Note that some variables in the csv raw data file might be in lowercase",
            ifelse(SD == ";", "          ;", ""), enter,
            "* To ensure that the dataset contains only variable names in uppercase",
            ifelse(SD == ";", "         ;", ""), enter, enter,
            "foreach var of varlist _all {",
            ifelse(SD == ";", " ;", ""), enter,
            "    local newname = upper(\"`var'\")",
            ifelse(SD == ";", " ;", ""), enter,
            "    cap rename \`var\' \`newname\'",
            ifelse(SD == ";", " ;", ""), enter,
            "}",
            ifelse(SD == ";", " ;", ""), enter, enter, enter, sep="")
        
        
        if (any(unlist(stringvars))) {
            cat("* Recode string variables which have labels, to numeric variables", ifelse(SD == ";", " ;", ""), enter, enter, sep="")
            for (i in seq(length(stringvars))) {
                if (stringvars[[i]]) {
                    cat(paste("encode (", toupper(names(stringvars)[i]), "), generate (", toupper(names(stringvars)[i]), "_NUM)", ifelse(SD == ";", " ;", ""), enter, sep=""))
                    vals <- seq(length(intrnlbls2$vallab[[names(stringvars)[i]]]))
                        ## sort() is necessary because Stata automatically transforms
                        # character to numeric using the labels in alfabetical ascending order (TO VERIFY THAT!!)
                    names(vals) <- names(sort(intrnlbls2$vallab[[names(stringvars)[i]]]))
                    intrnlbls2$vallab[[names(stringvars)[i]]] <- vals
                    newname <- paste(toupper(names(stringvars)[i]), "NUM", sep="_")
                    uniqueList2 <- lapply(uniqueList2, function(x) {
                        x[x == names(stringvars)[i]] <- newname
                        return(x)
                    })
                    names(intrnlbls2$vallab)[i] <- newname
                    names(intrnlbls2$varlab)[names(intrnlbls2$varlab) == names(stringvars)[i]] <- newname
                }
            }
            cat(enter, enter, sep="")
        }
        
        
        maxchars <- max(nchar(names(intrnlbls2$vallab)))
        
        cat("* Definition of variable labels", ifelse(SD == ";", " ;", ""), enter, enter, sep="")
        
        for (n in names(intrnlbls2$varlab)) {
            cat(paste("label variable ", toupper(n), paste(rep(" ", maxchars - nchar(n)), collapse=""), " \"", intrnlbls2$varlab[[n]][1], "\"", ifelse(SD == ";", " ;", ""), enter, sep=""))
        }
        cat(enter, enter, sep="")
        
        cat("* Definition of category labels", ifelse(SD == ";", " ;", ""), enter, enter, sep="")
        
        for (i in seq(length(uniqueList2))) {
            n <- uniqueList2[[i]][1]
            
            headerlabel <- paste("label define LABELS_GROUP_", i, sep="")
            
            if (SD == ";") {
                cat(headerlabel, enter,
                    paste(paste(paste(intrnlbls2$vallab[[n]], " \"", names(intrnlbls2$vallab[[n]]), "\"", sep=""),
                                collapse=enter), #paste("\n", paste(rep(" ", nchar(headerlabel)), collapse=""), sep="")),
                          " ;", enter, enter, sep=""),
                    sep="")
            }
            else if (SD == "///") {
                cat(headerlabel, " ///", enter,
                    paste(paste(paste(intrnlbls2$vallab[[n]], " \"", names(intrnlbls2$vallab[[n]]), "\"", sep=""),
                                collapse=" ///", enter), enter, enter, sep=""), sep="")
            }
            else if (SD == "/*") {
                cat(headerlabel, " /*", enter,
                    paste("*/ ", paste(paste(intrnlbls2$vallab[[n]], " \"", names(intrnlbls2$vallab[[n]]), "\"", sep=""),
                                collapse=paste(" /*", enter, "*/ ", sep="")),
                          enter, enter, sep=""),
                    sep="")
            }
            else {
                cat(paste("label define LABELS_GROUP_", i, " ", sep=""),
                    paste(paste(paste(intrnlbls2$vallab[[n]], " \"", names(intrnlbls2$vallab[[n]]), "\"", sep=""),
                                collapse=" "),
                          enter, sep=""),
                    sep="")
            }
        }
        
        cat(enter, enter, sep="")
        cat("* Attachment of category labels to variables", ifelse(SD == ";", " ;", ""), enter, enter, sep="")
        
        for (i in seq(length(uniqueList2))) {
            n <- uniqueList2[[i]][1]
            for (j in uniqueList2[[i]]) {
                cat(paste("label values ", toupper(j), paste(rep(" ", maxchars - nchar(j)), collapse=""), " LABELS_GROUP_", i, ifelse(SD == ";", " ;", ""), enter, sep=""))
            }
        }
        
        cat("\ncompress", ifelse(SD == ";", " ;", ""), enter,
            "save \"`statapath'\", replace", ifelse(SD == ";", " ;", ""), enter, enter,
            "log close", ifelse(SD == ";", " ;", ""), enter,
            "set more on", ifelse(SD == ";", "\n#delimit cr", ""), enter, sep="")
        
        sink()
        setwd(currentdir)
    }
    
    
    if (type == "SAS" | type == "all") {
        intrnlbls2 <- intrnlbls
        if (!file.exists("Setup files")) {
            dir.create("Setup files")
        }
        
        if (!file.exists(file.path("Setup files", "SAS"))) {
            dir.create(file.path("Setup files", "SAS"))
        }
        
        currentdir <- getwd()
        setwd(file.path("Setup files", "SAS"))
        sink(ifelse(length(grep("\\.sas", outfile)) > 0, outfile, paste(outfile, ".sas", sep="")))
        
        cat("* ------------------------------------------------------------------------------ ;", enter, enter,
            "* --- CONFIGURATION SECTION - START ---                                          ;", enter, enter, enter, sep="")                                            

        if (formats) {
            cat("* The following command should contain the complete path and                     ;", enter,
                "* name of the .csv file to be read (e.g. \"C:/CIS2008/Data/ALL.csv\")              ;", enter,
                "* Change CSV_DATA_PATH to your filename, below                                   ;", enter, enter,
                "FILENAME csvpath \"CSV_DATA_PATH\";", enter, enter, enter, sep="")
                      
           # cat("* It is assumed the data file was created under Windows (end of line is CRLF);", enter,
           #     "* If the csv file was created under Unix,  change eol=LF;", enter,
           #     "* If the csv file was created under MacOS, change eol=CR below;", enter, enter,
           #     "%LET eol=CRLF;", enter, enter, enter, sep="")
        }
        
        cat("* The following command should contain the complete path of the                  ;", enter,
            "* directory where the setup file will be saved (e.g. \"C:/CIS2008/Data\")          ;", enter,
            "* Change SAS_DATA_FOLDER to your directory name, below                           ;", enter, enter,
            "LIBNAME dirout \"SAS_DATA_FOLDER\";", enter, enter, enter, sep="")
                  
        cat("* The following command should contain the name of the output SAS file only      ;", enter,
            "* (without quotes, and without the .sas7bdat extension)                          ;", enter,
            "* Change SAS_FILE_NAME to your output file name, below                           ;", enter, enter,
            "%LET sasfile=SAS_FILE_NAME;", enter, enter, enter,
            "* --- CONFIGURATION SECTION -  END ---                                           ;", enter, enter,
            "* ------------------------------------------------------------------------------ ;", enter, enter, enter, enter,
            "* There should be nothing to change below this line;", enter,
            "* ------------------------------------------------------------------------------ ;", enter, enter, enter, enter, sep="")
        
        if (formats) {
            cat("* --- Read the raw data file ---                                                 ;", enter, enter,
                "DATA sasimport;", enter, enter,
                "INFILE csvpath", enter,
                "       DLM=", ifelse(delimiter == "\t", "'09'X", paste("\"", delimiter, "\"", sep="")), enter,
                "       FIRSTOBS=2", enter,
                # "       TERMSTR=&eol", enter,
                "       DSD", enter,
                "       TRUNCOVER", enter,
                "       LRECL=512", enter,
                "       ;", enter, enter,
                
                "INPUT  ", sasformats[1], toupper(csvnames[1]), enter, sep="")
            
            for (i in seq(2, length(csvnames))) {
                # cat("      ", sasformats[i], toupper(csvnames[i]), enter, sep="")
                cat("      ", toupper(csvnames[i]), sasformats[i], enter)
            }
            cat("       ;", enter)
            cat("RUN;", enter, enter,
                "* ------------------------------------------------------------------------------ ;", enter, enter, enter, sep="")
        }
        
        if (any(unlist(stringvars))) {
            cat("* --- Recode string variables which have labels, to numeric variables ---        ;", enter, enter,
                "DATA sasimport;", enter, enter,
                "    SET sasimport;", enter, enter, sep="")
            
            for (i in names(stringvars)) {
                if (stringvars[[i]]) {
                    x <- intrnlbls2$vallab[[i]]
                    vals <- seq(length(x))
                    names(vals) <- names(x)
                    for (j in seq(length(vals))) {
                        cat("    IF (", i, " = '", x[j], "') THEN ", i, " = ", vals[j], ";", enter, sep="")
                    }
                    
                    
                    cat(enter, "    TEMPVAR = input(", i, ", best12.);", enter,
                               "    DROP ", i, ";", enter,
                               "    RENAME TEMPVAR = ", i, ";", enter, enter, sep="")
                    
                    intrnlbls2$vallab[[i]] <- vals
                }
            }
            
            cat("RUN;", enter, enter,
                "* ------------------------------------------------------------------------------ ;", enter, enter, enter,
                
                "* --- Reorder the variables in the original positions ---                        ;", enter, enter,
                "DATA sasimport;", enter, enter,
                "    RETAIN ", gsub(",", "", splitrows(names(intrnlbls2$varlab), enter, 70, "           ")), ";", enter, enter,
                "    SET sasimport;", enter, enter,
                "RUN;", enter, enter,
                "* ------------------------------------------------------------------------------ ;", enter, enter, enter, sep="")
        }
        
        cat("* --- Add variable labels ---                                                    ;", enter, enter,
            "DATA sasimport;", enter, enter,
            "    SET sasimport;", enter, enter, sep="")
        
        for (i in seq(length(varnames))) {
            cat("    LABEL ", toupper(varnames[i]), paste(rep(" ", maxchars - nchar(varnames[i]) + 1), collapse=""), "=", " \"", intrnlbls2$varlab[[i]][1], "\";", enter, sep="")
        }
        
        cat(enter, "RUN;", enter, enter,
            "* ------------------------------------------------------------------------------ ;", enter, enter, enter, sep="")
        
        cat("* --- Create value labels groups ---                                             ;", enter, enter,
            "PROC FORMAT;", enter, enter, sep="")
        
        for (i in seq(length(uniqueList))) {
            n <- uniqueList[[i]][1]
            cat(paste("VALUE LABELS_", i, "_GROUP", enter, sep=""),
                paste(paste(paste("      ", intrnlbls2$vallab[[n]], "=\"", names(intrnlbls2$vallab[[n]]), "\"", sep=""),
                            collapse=enter), #paste("\n", paste(rep(" ", nchar(headerlabel)), collapse=""), sep="")),
                      ";", enter, enter, sep=""),
                sep="")
        }
		  
        cat("RUN;", enter, enter,
            "* ------------------------------------------------------------------------------ ;", enter, enter, enter, sep="")
        
        cat("* --- Format variables with value labels ---                                     ;", enter, enter,
            "DATA sasimport;", enter, enter,
            "    SET sasimport;", enter, enter, "    FORMAT", enter, sep="")
        
        for (i in seq(length(uniqueList))) {
            n <- uniqueList[[i]][1]
            for (j in uniqueList[[i]]) {
                cat("    ", toupper(j), paste(rep(" ", maxchars - nchar(j)), collapse=""), " LABELS_", i, "_GROUP", ".", enter, sep="")
            }
        }
          
        cat("    ;", enter, enter, "RUN;", enter, enter,
            "* ------------------------------------------------------------------------------ ;", enter, enter, enter, sep="")
        
                      
        cat("* --- Save data to a sas type file ---                                           ;", enter, enter,
            "DATA dirout.&sasfile;", enter, enter,
            "    SET sasimport;", enter, enter,
            "RUN;", enter, sep="")
        
        sink()
        setwd(currentdir)
    }
    
    
    
    if (type == "R" | type == "all") {
        intrnlbls2 <- intrnlbls
        printMISSING <- FALSE
        
        if (!file.exists("Setup files")) {
            dir.create("Setup files")
        }
        
        if (!file.exists(file.path("Setup files", "R"))) {
            dir.create(file.path("Setup files", "R"))
        }
        
        currentdir <- getwd()
        setwd(file.path("Setup files", "R"))
        sink(ifelse(length(grep("\\.R", outfile)) > 0, outfile, paste(outfile, ".R", sep="")))
        
        cat("# ------------------------------------------------------------------------------", enter, enter,
            "# --- CONFIGURATION SECTION - START ---", enter, enter, enter, sep="")

        if (formats) {
            cat("# The following command should contain the complete path and", enter,
                "# name of the .csv file to be read (e.g. \"C:/CIS 2008/Data/ALL.csv\")", enter,
                "# Change CSV_DATA_PATH to your filename, below:", enter, enter,
                "csvpath <- \"CSV_DATA_PATH\"", enter, enter, enter, sep="")
        }
        
        cat("# The following command should contain the complete path and", enter,
            "# name of the .Rdata file to be saved (e.g. \"C:/CIS 2008/Data/ALL.Rdata\")", enter,
            "# Change RDATA_PATH to your filename, below:", enter, enter,
            "rdatapath <- \"RDATA_PATH\"", enter, enter, enter,
            # "# The following command should contain the name of the unique ID variable", enter,
            # "# in the .csv file (to identify missing observations in different variables).", enter,
            # "# Change UNIQUEID to your unique ID variable name, below:", enter, enter,
            # "uniqueid <- \"UNIQUEID\"", enter, enter, enter,
            "# --- CONFIGURATION SECTION -  END  ---", enter, enter,
            "# ------------------------------------------------------------------------------", enter, enter, enter, enter,
            "# There should be nothing to change below this line", enter,                                                            
            "# ------------------------------------------------------------------------------", enter, enter, enter, enter, sep="")
                  
        if (formats) {
            cat("# --- Read the raw data ---", enter, enter,
                "rdatafile <- read.csv(csvpath, sep = \"", ifelse(delimiter == "\t", "\\t", delimiter), "\")",
                enter, enter, "names(rdatafile) <- toupper(names(rdatafile))    # all variable names to upper case",
                enter, enter, enter,
                "# ------------------------------------------------------------------------------",
                enter, enter, enter, enter, sep="")
        }
        else {
            cat("# \"rdatafile\" should be an R data.frame (usually read from a .csv file)\n\n")
        }
        
        
        if (any(unlist(stringvars))) {
            cat("# --- Recode string variables which have labels, to numeric variables ---", enter, enter, sep="")
            stringvars <- stringvars[which(unlist(stringvars))]
            for (i in names(stringvars)) {
                if (stringvars[[i]]) {
                    x <- intrnlbls2$vallab[[i]]
                    vals <- seq(length(x))
                    names(vals) <- names(x)
                    cat(i, " <- rdatafile[ , \"", i, "\"]", enter,
                        "tempvar <- rep(NA, length(", i, "))", enter, enter, sep="")
                    for (j in seq(length(vals))) {
                        cat("tempvar[", i, " == \"", x[j], "\"] <- ", vals[j], enter, sep="")
                    }
                    
                    cat(enter, "rdatafile[ , \"", i, "\"] <- tempvar", enter, enter, sep="")
                    
                    intrnlbls2$vallab[[i]] <- vals
                }
            }
            
            cat(enter,
                "# ------------------------------------------------------------------------------",
                enter, enter, enter, enter, sep="")
        }
        
        cat("# --- Set the variable labels attribute --- ", enter, enter,
            "attr(rdatafile, \"variable labels\") <- list(", enter, sep="")
            
        for (i in varnames) {
            cat("\"", toupper(i), "\"",
            paste(rep(" ", maxchars - nchar(i)), collapse=""),
            " = \"", intrnlbls2$varlab[[i]][1], "\"",
            ifelse(i == varnames[length(varnames)], ")", ","),
            enter, sep="")
        }
        
        uList <- unlist(uniqueList)
        
        cat(enter, enter,
            "# ------------------------------------------------------------------------------",
            enter, enter, enter, enter,
            "# --- Set the value labels attribute --- ", enter, enter,
            "attr(rdatafile, \"value labels\") <- vector(mode=\"list\", length=", length(uList), ")", enter, enter, sep="")
            # "names(attr(rdatafile, \"value labels\")) <- c(", enter, 
            # splitrows(paste("\"", uList, "\"", sep=""), enter, 80), enter,
            # ")", enter, enter, sep="")
        
        contor <- 1
        
        for (i in seq(length(uniqueList))) {
            n <- uniqueList[[i]][1]
            
            listelements <- ifelse(length(uniqueList[[i]]) == 1, contor, paste(contor, ":", contor + length(uniqueList[[i]]) - 1, sep=""))
            
            cat("attr(rdatafile, \"value labels\")[", listelements, "] <- ", 
                ifelse(length(uniqueList[[i]]) == 1, "list(c(", "rep(list(c("), enter,
                paste(paste("    \"", names(intrnlbls2$vallab[[n]]), "\" =", sep=""),
                      intrnlbls2$vallab[[n]],
                      collapse=paste(",", enter, sep="")),
                sep="")
            
            cat(enter, ifelse(length(uniqueList[[i]]) == 1, "))", 
                paste(")), ", length(uniqueList[[i]]), ")", sep="")), enter, enter,
                "names(attr(rdatafile, \"value labels\"))[", listelements, "] <- ", sep="")
            if (length(uniqueList[[i]]) == 1) {
                cat("\"", uniqueList[[i]], "\"", enter, enter, sep="")
            }
            else {
                cat("c(", enter, "    ",
                    splitrows(paste("\"", uniqueList[[i]], "\"", sep=""), enter, 80, spacerep="    "), enter,
                    ")", enter, enter, sep="")
            }
            
            contor <- contor + length(uniqueList[[i]])
        }
        
        
        cat(enter,
            "# ------------------------------------------------------------------------------",
            enter, enter, enter, enter, sep="")
        
        
        printMISSING <- FALSE
        if (!missing(miss) & uniqueid != "" & is.data.frame(csv)) {
            
            if (is.numeric(miss)) {
                missvars <- lapply(intrnlbls2$vallab, function(x) !is.na(match(miss, x)))
                withmiss <- as.vector(unlist(lapply(missvars, any)))
                missvals <- lapply(missvars[withmiss], function(x) miss[x])
                names(missvals) <- names(intrnlbls2$vallab)[withmiss]
            }
            else {
                missvars <- lapply(intrnlbls2$vallab, function(x) !is.na(match(names(x), miss)))
                withmiss <- as.vector(unlist(lapply(missvars, any)))
                missvals <- lapply(which(withmiss), function(x) as.vector(intrnlbls2$vallab[[x]][missvars[[x]]]))
                names(missvals) <- names(intrnlbls2$vallab)[withmiss]
            }
            
            
            cat("# --- Set the missing values attribute --- ", enter, enter,
                "attr(rdatafile, \"unique id\") <- \"", uniqueid, "\"", enter, enter, sep="")
                
            for (i in seq(length(missvals))) {
                values <- intrnlbls2$vallab[[names(missvals)[i]]]
                cat("# ", names(missvals)[i], enter, sep="")
                
                for (j in seq(length(missvals[[i]]))) {
                    
                    # for the moment, the length of missvals[[i]][j] is always equal to 1,
                    # but in the future metadata information might allocate ranges for some missing types
                    # therefore missvals[[i]] might be a list
                    
                    testvals <- ifelse(length(missvals[[i]][j]) == 1,
                                       missvals[[i]][j],
                                       paste("c(", paste(missvals[[i]][j], collapse = ", "), ")", sep=""))
                    
                    cat("attr(rdatafile, \"missing types\")$", names(missvals)[i],
                        "[[\"", names(values)[values == missvals[[i]][j]], "\"]] <- list(", enter,
                        "values = ", testvals, ",", enter,
                        "cases = rdatafile$", uniqueid, "[rdatafile$", names(missvals)[i], 
                        ifelse(length(missvals[[i]][j]) == 1, " == ", " %in% "), testvals,  
                        "]", enter, ")", enter, enter, sep="")
                }
                
                # here, though, it is likely to have multiple missing values
                # in the future, missvals[[i]] might be a list to accommodate for ranges
                
                cat("rdatafile$", names(missvals)[i], "[rdatafile$", names(missvals)[i],
                    ifelse(length(unlist(missvals[[i]])) == 1, paste(" == ", unlist(missvals[[i]]), sep=""),
                           paste(" %in% c(", paste(unlist(missvals[[i]]), collapse = ", "), ")", sep="")), 
                    "] <- NA", enter, enter, sep="")
            }
            
            cat(enter, enter,
            "# ------------------------------------------------------------------------------",
            enter, enter, enter, enter, sep="")
        }
        else if (!missing(miss) & (uniqueid == "" | is.data.frame(csv))) {
            printMISSING <- TRUE
        }
        
        cat("# --- Save the R data file --- ", enter, enter,
            "rfilename <- unlist(strsplit(basename(rdatapath), split=\"\\\\.\"))[1]", enter,
            "rdatapath <- file.path(dirname(rdatapath), paste(rfilename, \".Rdata\", sep=\"\"))", enter,
            "assign(rfilename, rdatafile)", enter,
            "eval(parse(text = paste(\"save(\", rfilename, \", file=rdatapath)\", sep=\"\")))",
            enter, enter, enter,
            "# ------------------------------------------------------------------------------",
            enter, enter, enter, enter,
            "# --- Clean up the working space --- ", enter, enter,
            "rm(rfilename, rdatafile, csvpath, rdatapath, tempvar, vallist", sep="")
        
        if (any(unlist(stringvars))) {
            cat(", ", paste(names(stringvars)[unlist(stringvars)], collapse=", "), sep="")
        }
        
        cat(")", enter, enter, enter,
            "# ------------------------------------------------------------------------------",
            enter, enter, enter, enter, sep="")
        
        
        
        # finish writing and close the .R file
        sink()
        
        if (printMISSING) {
            # this would be printed on the screen
            cat("The \"csv\" and \"uniqueid\" arguments are both mandatory to produce R missing values commands.\n\n")
        }
                
        setwd(currentdir)
        
    }
}

