##******************************************************************************
## GUI for the Renext package.
##
## Coding author:  Yves Deville and IRSN.
##
## Design and specification: Lise Bardet and Claire-Marie Duluc (IRSN) 
##
## Please setup your working directory to the place where csv files
## will be found by default
##
##******************************************************************************

## Note YD 2011-04-14
##
## Les commandes "enabled <- " produisent des erreurs et empechent
## l'utilisation du package!!!
## Cela semble lie a la commande "state" de Tcl qui n'est pas
## trouvee.
## 
## Le pb se produit ave R-2.12, R-2.13, mais pas en version R-2.11


## .onLoad <- function (libname, pkgname) {
##  cat("appel de .onLoad\n")
##  dataList <<- makeDemos()
##  fitsList <<- list("<New>" = list(coef = NA, pred = NULL))
##  resList <<- list()
##}

RenextGUI <- function(guiToolkit = c("tcltk", "RGtk2", "WWW"),
                      useGlobalenv = FALSE,
                      debug = FALSE)  {

  ## Create a new environment where most things will be stored?
  
    if (useGlobalenv) {
        re <- globalenv()
    } else {
        re <- environment()
    }
    
    DEBUG <- debug
    strW <- paste(rep(" ", 40), collapse = "")
    
    ## require(gWidgets)
    ## require(Renext)
    
    if (DEBUG) cat("initialisation tasks...")
    
    guiToolkit <- match.arg(guiToolkit)
    options(guiToolkit = guiToolkit)
    
    dataOK <- FALSE
    main <- ""
    
    ## define
    odtModel <- file.path(system.file("reporting", package = "RenextGUI"),
                          "ExIn.odt")
    
    reportingDir <- getwd() 
    
    ## List for storing temporary objects or values 
    re$RLprov <- list()
    
    
    re$RLprov$.decs <- c(".", ",")
    names(re$RLprov$.decs) <- c("\".\"(dec. point)", "\",\" (comma)")
    
    ## column separators for vsc files
    re$RLprov$.seps <- c(";", "\t", " ", ",")
    names(re$RLprov$.seps) <- c("\";\" (semi-column)", "\"\\t\" (tab)",
                                "\" \"(space)", "\",\" (comma)")
    
    ## list to handle widgets displayed in the results noteBook
    ## items are the widgets paths
    re$.nB <- list()
    
    ##============================================================================
    ## default values displayed when no project or no fit is selected.
    ##
    ## Note that in ".defaultFit", there are no 'rMax1' nor 'rMax2' items
    ## since these widgets do not exist on the 'fit' item
    ##============================================================================
    
    re$RLprov$.defaultProj <-
        list(file = " ",
             varName = " ",
             main = " ",
             effDuration = "",
             rMax1 = " ",
             rMax1Dur = " ",
             rMax2 = " ",
             rMax2Dur = " ")
    
    re$RLprov$.defaultFit <-
        list(distribution = "exponential",
             useMax1 = "no",
             rMax1Lab = "                                   ",
             useMax2 = "no",
             rMax2Lab = "                                   ",
             threshold = 0,
             effDuration = "",
             confLev = "95;70",
             returnPeriods = "10;20;50;100;200;500;1000")
    
    re$dataList <- makeDemos()
    re$fitsList <- list()
    re$resList <- list()
    
    ##****************************************************************************
    ## miscellanea functions for use in RenextGui devoted to 'projects'
    ##
    ##
    ## Note that these function make access to 'global' objects which are 
    ## NOT in the scope of the function.
    ##
    ## These are widgets of lists containing data.
    ## 
    ##****************************************************************************
    
    ##============================================================================
    ## keep only project names with valid data
    ##============================================================================
    
    ## cleanProjects <- function(h, ...) {
    ##   
    ## }
    
    ##============================================================================
    ## Plot of a dataset
    ##
    ## The handler is unimportant here.
    ##
    ##============================================================================
    
    plotData <- function(h, ...) {
        
        if (DEBUG) cat("[ Entering 'plotData'\n")
        
        project <- svalue(datW$project)
        
        if (DEBUG) cat("project = ", project, "\n")
        
        if (length(project) && !is.na(project)) {
            
            item <- re$dataList[[project]]
            
            if (DEBUG) {
                cat("varName = ", item$varName, "\n")
                cat("names(item = ", names(item))
            }
            
            if (item$hasDate) {
                plot(x = item$data$date,
                     y = item$data$x,
                     col = "orangered", bg = "orange",
                     type = "h", pch = 21,
                     ylab = item$varName,
                     xlab = "date",
                     main = svalue(datW[["main"]]))
            } else {
                plot(x = item$data$x,
                     col = "orangered", bg = "orange",
                     type = "h", pch = 21,
                     ylab = item$varName,
                     xlab = "obs.",
                     main = svalue(datW[["main"]]))
            }
            
        }
        if (DEBUG) {
            cat("] Exiting 'plotData'\n")
        }
    }
    
    ##=====================================================================
    ## Change the definition of a dataset
    ##
    ##=====================================================================
    
    updateData <- function(h, ...) {

        if (DEBUG) {
            cat("[ Entering 'updateData'\n")
            cat("names(h$action) = ", names(h$action), "\n")
        }
        
        svalue(datW$sb) <- ""
        
        if (length(h$obj[])) {
            project <- svalue(h$obj)
        } else {
            project <- NA
        }
        
        if (DEBUG) cat("project :   ", project, "\n")
        
        if ( length(project) || !is.na(project) ) {
            
            item <- re$dataList[[project]]
            
            if (DEBUG) {
                cat("names(item)\n")
                print(names(item))
            }
            
            ## CANCELED: DOES NOT WORK WITH tcltjk
            ## 
            ## for (name in c("file", "main", "varName", "effDuration",
            ##               "rMax1", "rMax1Dur", "rMax2", "rMax2Dur")) {
            ##  if ( !is.null(datW[[name]]) ) { 
            ##    if ( length(item[[name]]) && !is.na(item[[name]]) ) { 
            ##      svalue(datW[[name]]) <- item[[name]]
            ##    } else {
            ##      svalue(datW[[name]]) <- ""
            ##    } 
            ##  }    
            ## }
            svalue(datW[["file"]]) <- item[["path"]]
            svalue(datW[["main"]]) <- item[["main"]]
            svalue(datW[["varName"]]) <- item[["varName"]]
            svalue(datW[["effDuration"]]) <- item[["effDuration"]]
            svalue(datW[["rMax1"]]) <- item[["rMax1"]]
            svalue(datW[["rMax1Dur"]]) <- item[["rMax1Dur"]]
            svalue(datW[["rMax2"]]) <- item[["rMax2"]]
            svalue(datW[["rMax2Dur"]]) <- item[["rMax2Dur"]]
            
            ## MODIFS 2011-12-13 added the test on 'isLocked'
            if ( is.null(item$isLocked) || !item$isLocked ) {
                
                ## DOES NOT WORK AS EXPECTED
                for ( name in names(re$RLprov$.defaultProj) ) {
                    if (DEBUG) cat("name = ", name, "is ENabled\n")
                    enabled(datW[[name]]) <- TRUE
                }
                enabled(datW[["OK"]]) <- TRUE
                svalue(datW$sb) <- "edit project if needed"
                
            } else {
                
                for ( name in names(re$RLprov$.defaultProj) ) {
                    if (DEBUG) cat("name = ", name, "is DISabled\n")
                    enabled(datW[[name]]) <- FALSE
                }
                enabled(datW[["OK"]]) <- FALSE
                
                svalue(datW$sb) <- "project can not be edited"
            }
            
        } else {
            
            if (TRUE) {
                ## MODIFS 2011-10-14 
                ## reset to default values
                for ( name in names(re$RLprov$.defaultProj) ) {
                    svalue(datW[[name]]) <- re$RLprov$.defaultProj[[name]]
                    enabled(datW[[name]]) <- FALSE
                }
                
                enabled(datW[["OK"]]) <- FALSE
                
                svalue(h$action$sb) <- "no project selected"
                
                ## END MODIFS
            } else {
                
                ## for (name in c("file", "main", "varName")) {
                ##  svalue(datW[[name]]) <- ""
                ## }
                svalue(datW[["file"]]) <- ""
                svalue(datW[["main"]]) <- ""
                svalue(datW[["varName"]]) <- ""
                
            }
            
        }
        
        if (DEBUG) {
            cat("] Exiting 'updateData'\n")
        }
        
    }
    
    ##=======================================================================
    ## Modify 'dataList' from the widgets svalue 
    ##=======================================================================
    
    bindData <- function(h, ...) {
        
        if (DEBUG) cat("[ Entering 'bindData'\n")
        
        project <- svalue(datW$project)
        if (DEBUG) cat("project = ", project, "\n")
        
        if ( !is.na(project) && (length(project) > 0) ) {
            
            name <- h$action
            if (DEBUG) cat("name = ", name, "\n")
            ## 2011-12-14 added the gsub
            re$dataList[[project]][[name]] <- gsub(" ", "", svalue(datW[[name]]))
            if (DEBUG) cat("***",  re$dataList[[project]][[name]], "***\n")
        } 
        
        if (DEBUG) cat("] Exiting 'bindData'\n")
        
    }

    ##============================================================================
    ## Modify 'dataList' from the widgets svalue 
    ##============================================================================
    
    bindAllData <- function(h, ...) {
        
        
        if (DEBUG) cat("[ Entering 'bindAllData'\n")
        
        project <- svalue(datW$project)
        name <- h$action
        
        ## 2011-12-14 added the gsubs
        re$dataList[[project]][["main"]] <-
            gsub(" ", "", svalue(datW[["main"]]))
        re$dataList[[project]][["effDuration"]] <-
            gsub(" ", "",  svalue(datW[["effDuration"]]))
        re$dataList[[project]][["rMax1"]] <-
            gsub(" ", "", svalue(datW[["rMax1"]]))
        re$dataList[[project]][["rMax1Dur"]] <-
            gsub(" ", "", svalue(datW[["rMax1Dur"]]))
        re$dataList[[project]][["rMax2"]] <-
            gsub(" ", "", svalue(datW[["rMax2"]]))
        re$dataList[[project]][["rMax2Dur"]] <-
            gsub(" ", "", svalue(datW[["rMax2Dur"]]))
        
        if (DEBUG) cat("] Exiting 'bindAllData'\n")
        
    }
    
    ##==========================================================================
    ## NEW IN VERSION 0.5-1
    ##
    ## Try to read csv file and find suitable variables:
    ##
    ## o one date variable named 'date', and
    ## o one numeric variable.
    ## 
    ## The parameters can be edited by user and will be stored into the
    ## dataList item created if any. 
    ##  
    ##==========================================================================

    tryReadCsv <- function(h, ...) {
        
        if (DEBUG) {
            cat("entering 'tryReadCsv'\n")
            cat("source = ", h$action$source, "\n")
        }
        
        svalue(h$action$sb) <- "changing file parameters"
        
        fileName <- svalue(h$action$file)
        
        if (DEBUG) cat(sprintf("fileName = %s\n", fileName))
        
        re$RLprov$projectPar$path <- fileName
        
        ##---------------------------------------------------------------------
        ## if file is changed (or at least refreshes), reinitialize
        ## the proj. name
        ## --------------------------------------------------------------------
    
        if (!is.null(h$action$source) && (h$action$source == "file") ) {
            
            if ( !file.exists(fileName) ) {
                svalue(h$action$sb) <- "select an existing file"
                return() 
            }
            
            ## the tail part of the path
            fileNameShort <- rev(unlist(strsplit(as.character(fileName),
                                                 split = .Platform$file.sep)))[1]
            
            ## guess name(s)
            fs <- unlist(strsplit(as.character(fileNameShort), split = "\\."))
            if ((nc <- length(fs)) > 1) fs <- paste(fs[1:(nc-1)], collapse = ".")
            
            projName <- guessNames(cand = fs, existing = names(re$dataList))
            if (DEBUG) cat("Setting project name to ", projName, "\n")
            svalue(h$action$project) <- projName
            
        }
        
        ##======================================================================
        ## If we need read the file... that is if we do not simply
        ## change the variables columns or names
        ##======================================================================
        
        if (!is.null(h$action$source) &&
            !(h$action$source %in% c("dateCol", "varCol", "varName")) ) {
            
            sep <- re$RLprov$.seps[as.integer(svalue(h$action$sep, index = TRUE))]
            if (DEBUG) cat("sep = ", sep, "\n")
            
            header <- svalue(h$action$header)
            skip <- svalue(h$action$skip)
            if (length(skip) == 0L) {
                skip <-   svalue(h$action$skip) <- 0L
                
            }
            
            ## for later ???
            if (!is.null(h$action$dec)) {
                dec <- re$RLprov$.decs[as.integer(svalue(h$action$dec, index = TRUE))]
            } else {
                dec <- "."
            }
            
            if (DEBUG) {
                cat("sep = ", sep, "header = ", header, "skip = ",
                    skip, "dec = ", dec, "\n")
            }
            test <- try(readLines(fileName, n = 7L))
            
            if (is(test, "try-error")) {
                svalue(h$action$text) <- "<first lines>"
                alarm() 
                svalue(h$action$sb) <- geterrmessage()
            } else {
                nf <- count.fields(file = fileName, sep = sep)
                svalue(h$action$text) <- test
            }
            
            res <- try(read.table(file = fileName,
                                  header = header,
                                  skip = skip,
                                  sep = sep,
                                  stringsAsFactors = FALSE,
                                  dec = dec))
            
            if (is(res, "try-error")) {
                cat("ERROR\n")
                alarm()
                h$action$data <- NULL
                svalue(h$action$sb) <- geterrmessage()
                enabled(h$action$OK) <- FALSE
            } else {
                
                dat <- head(res)
                
                re$RLprov$projectPar$data <- res
                
                ##=============================================================
                ## Find out which columns are date and which are numeric
                ##=============================================================
                
                isDate <- rep(FALSE, ncol(dat))
                isNum <- rep(FALSE, ncol(dat))
                
                dateFormat <- rep("<none>", ncol(dat))
                ## for each column: can it be a date (POSIXct)?
                
                for (i in 1:ncol(dat)) {
                    
                    if (is.character(dat[[i]])) {
                        ## cat("i = ", i, "character\n")
                        
                        if (all(regexpr("[0-2][0-9][0-9][0-9]-[0-1][0-9]-[0-3][0-9]$",
                                        dat[[i]])>0) ) {
                            if (DEBUG) cat("i = ", i, "POSIXct\n")
                            res[[i]] <- as.POSIXct(res[[i]])
                            isDate[i] <- TRUE
                            dateFormat[i] <- "%Y-%m-%d" 
                        } else if (all(regexpr("[0-2][0-9][0-9][0-9]/[0-1][0-9]/[0-3][0-9]$",
                                               dat[[i]])>0) ) {
                            res[[i]] <- as.POSIXct(strptime(res[[i]], "%Y/%m/%d"))
                            isDate[i] <- TRUE
                            dateFormat[i] <- "%Y/%m/%d" 
                        } else if (all(regexpr("[0-3][0-9]/[0-1][0-9]/[0-2][0-9][0-9][0-9]$",
                                               dat[[i]])>0) ) {
                            res[[i]] <- as.POSIXct(strptime(res[[i]], "%d/%m/%Y"))
                            isDate[i] <- TRUE
                            dateFormat[i] <- "%d/%m/%Y"
                        } 
                        
                    } else if (is.numeric(dat[[i]])) {
                        if (DEBUG) cat("i = ", i, "numeric\n")
                        isNum[i] <- TRUE
                    }
                    
                    h$action$data <- res
                    
                }
                
                ## stor nums to the temporary list within 'action'
                re$RLprov$projectPar$nCol <- ncol(res)
                re$RLprov$projectPar$indDate <- (1L:ncol(res))[isDate]
                re$RLprov$projectPar$indNum <- (1L:ncol(res))[isNum]
                
                
                ##=============================================================
                ## Using the previous results, populate the selection
                ## lists for the date and for the numeric variable.
                ## ============================================================
      
                if (any(isNum)) {
                    
                    if (DEBUG) cat("at least one 'num' column found\n")
                    
                    nms <- paste("(", paste("col", 1:ncol(res), sep ="_"), ")", sep = "")
                    ## svalue(h$action$varCol) <- ""
                    h$action$dateCol[] <- paste(colnames(res)[isDate], nms[isDate],
                                                sep = "_")
                    h$action$varCol[] <- paste(colnames(res)[isNum], nms[isNum],
                                               sep = "_")
                    
                    svalue(h$action$sb) <- "choose variable(s) and rename if needed"
                    
                    ## if ther is at least on date column, enable choice
                    if (any(isDate)) {
                        
                        if (DEBUG) cat("at least one 'date' column found\n")
                        
                        svalue(h$action$dateFormat, index = TRUE) <- dateFormat[isDate][1L]
                        svalue(h$action$dateCol, index = TRUE) <- 1L
                        re$RLprov$projectPar$hasDate <- TRUE
                        
                    } else {
                        
                        if (DEBUG) cat("no 'date' col found\n")
                        
                        h$action$dateCol[] <- ""
                        svalue(h$action$dateCol) <- ""
                        re$RLprov$projectPar$hasDate <- FALSE
                        
                    }
                    
                    svalue(h$action$varCol, index = TRUE) <- 1L
                    svalue(h$action$varName) <- colnames(res)[isNum][1L]
                    
                    if (DEBUG) {
                        cat("enable \"OK\" button\n")
                    }
                    enabled(h$action$OK) <- TRUE
                    svalue(h$action$sb) <- sprintf("data OK with %d rows", nrow(res))
                    
                } else {
                    
                    if ( options()$guiToolkit == "RGtk2" ) {
                        h$action$dateCol[] <- character(0)
                        h$action$varCol[] <- character(0)
                        svalue(h$action$dateCol, index = TRUE) <- 0L
                        svalue(h$action$varCol, index = TRUE) <- 0L
                    } else {
                        h$action$dateCol[] <- ""
                        h$action$varCol[] <- ""
                        svalue(h$action$dateCol, index = TRUE) <- 1L
                        svalue(h$action$varCol, index = TRUE) <- 1L
                    }
                    ## alarm()
                    svalue(h$action$sb) <-
                        "no numeric variable found: try different parameters"
                    svalue(h$action$dateFormat) <- "<none>"
                    svalue(h$action$varName) <- " "
                    enabled(h$action$OK) <- FALSE
                    
                    re$RLprov$projectPar$hasDate <- FALSE
                    ## RLprov$projectPar$csvPar$dateColNum <- NA
                    ## RLprov$projectPar$csvPar$varColNum <- NA
                    
                }
                
                if (DEBUG) {
                    cat("svalue dateCol", svalue(h$action$dateCol), "\n")
                    cat("svalue varCol", svalue(h$action$varCol), "\n")
                    ## print(h$action$data)
                }
                
            }
            ## cat("classes = \n")
            ## print(lapply(dat, function(x) class(x)[1]))
            
        } else if (h$action$source == "varCol") {
            
            mots <- unlist(strsplit(svalue(h$action$varCol), split = "_"))
            svalue(h$action$varName) <- mots[1]
            ## RLprov$projectPar$csvPar$varColNum <-
            ##  as.integer(substring(mots[3], first = 1L, last = nchar(mots[3])- 1L))
            
            
        }  else if (h$action$source == "dateCol") {
            
            mots <- unlist(strsplit(svalue(h$action$dateCol), split = "_"))
            svalue(h$action$dateName) <- mots[1]
            
        }
        
        re$RLprov$projectPar$varName <- svalue(h$action$varName)
        re$RLprov$projectPar$sep <- svalue(h$action$sep)
        re$RLprov$projectPar$header <- svalue(h$action$header)
        re$RLprov$projectPar$skip <- svalue(h$action$skip)
        re$RLprov$projectPar$dateFormat <- svalue(h$action$dateFormat)
        
        
        if (DEBUG) {
            cat("Parameters\n")
            str(re$RLprov$projectPar)
        }
        
    }
    
    ##==========================================================================
    ## Dialogs to validate a new project from a csv file
    ##
    ## Here the hander 'h' has
    ##
    ## - 'obj'      the dialog window,
    ## - 'action'   the list of widgets in the dialog window
    ##
    ## The values of 'Main', 'Path', are found from the widgets
    ## and copied in the 'dataList' object in the parent env(!)
    ## 
    ##============================================================================
    
    validNewCsvProject <- function(h, ...) {
        
        if (DEBUG) cat("entering 'validNewCsvProject'\n")
        ## widgets List
        wL <- h$action
        
        if (DEBUG) cat("names(wL)", names(wL), "\n")
        
        projectName <- svalue(wL$project)
        varName <- as.character(svalue(wL$varName))
        
        if ( is.na(varName) ) {
            alarm( )
            svalue(wL$sb) <- "Bad variable name"
            return()
        }
        
        if ( projectName %in% names(re$dataList) ) {
            alarm()
            svalue(wL$sb) <- "Project name already exists!!! Please change it."
        } else {
            
            if (DEBUG) cat("new project name...\n")
            
            colNum <- as.integer(svalue(wL$varCol, index = TRUE))
            colNum <-  re$RLprov$projectPar$indNum[colNum]
            
            if (re$RLprov$projectPar$hasDate) {
                
                colDate <- as.integer(svalue(wL$dateCol, index = TRUE))
                colDate <-  re$RLprov$projectPar$indDate[colDate]
                
                ## Modif 2011-12-13 
                dat <- data.frame(date = strptime(re$RLprov$projectPar$data[ , colDate],
                                      format = re$RLprov$projectPar$dateFormat),
                                  ##date = as.POSIXct(RLprov$projectPar$data[ , colDate]),
                                  x = re$RLprov$projectPar$data[ , colNum])
                
            } else {
                colDate <- NA
                dat <- data.frame(x = re$RLprov$projectPar$data[ , colNum])
            }
            
            newItem <-
                list(source = "csv file",
                     csvPar = list(nCol = re$RLprov$projectPar$nCol,
                         header = re$RLprov$projectPar$header,
                         sep = re$RLprov$.seps[as.integer(svalue(wL$sep, index = TRUE))],
                         skip = as.integer(svalue(wL$skip)),
                         dec = re$RLprov$.decs[as.integer(svalue(wL$dec, index = TRUE))]),
                     hasDate = re$RLprov$projectPar$hasDate,
                     colDate = colDate,
                     dateFormat =  re$RLprov$projectPar$dateFormat,
                     colNum = colNum,
                     path = svalue(wL$file),
                     info = list(shortLab = svalue(wL$project),
                         varName = svalue(wL$varName)),
                     data = dat,
                     varName = svalue(wL$varName),
                     main = projectName,
                     effDuration = " ",
                     thresholds = guessThreshold(dat$x),
                     rMax1 = " ",
                     rMax1Dur = " ",
                     rMax2 = " ",
                     rMax2Dur = " ")
            
            if (DEBUG) str(newItem)
            
            
            svalue(datW$sb) <- "New project added"
            L <- length(re$dataList)
            re$dataList[[projectName]] <- newItem
            
            if (L > 0) {
                datW$project[1:(L+1)] <- c(datW$project[1:L], projectName)
            } else {
                datW$project[] <- projectName
            }
            
            ## Former action removed
            ## if ( options()$guiToolkit == "RGtk2" ) {
            ##   svalue(datW$project) <- projectName
            ## }
        }
        
        dispose(wL[["window"]])
        
        if (DEBUG) cat("exiting 'validNewCsvProject'\n")
        
    }
    
    ##=========================================================================
    ## Create a new project by reading a csv file for the data part.
    ## 
    ##=========================================================================
    
    newCsvProject <- function(h, ...) { 
        
        
        wFile <- gwindow("csv File Import Wizard",
                         height = 500, width = 500, visible = FALSE)
        
        fiW <- list( )
        
        ##====================================================================
        ## defaut items
        ##====================================================================
    
        re$RLprov$projectPar <- list(path = NULL,
                                     data = NULL,
                                     varName = character(0),
                                     indDate = integer(0),
                                     indNum = integer(0),
                                     sep = 1L,
                                     header = TRUE,
                                     skip = 0L,
                                     dateFormat = "%Y-%m-%d",
                                     hasDate = FALSE)
        
        fiTbl <- glayout(container = wFile, label = "Csv file import",
                         cexpand = TRUE, fill = "both")
        ##===================================================================
        fiTbl[1, 1] <- glabel("OT data csv file", container = fiTbl)
        
        fiTbl[1, 2:6] <-
            ( fiW$file <- gfilebrowse("File_to_read",
                                      container = fiTbl,
                                      width = 50,
                                      quote = FALSE) )
        ##===================================================================
        fiTbl[2, 1] <- glabel("Project", container = fiTbl)
        fiTbl[2, 2:6] <- ( fiW$project <- gedit("                 ",
                                                container = fiTbl) )
        
        ##===================================================================
        fiTbl[3, 1:2] <- glabel("First lines in file", container = fiTbl)
        ##===================================================================
        fiTbl[4, 1:6] <- (fiW$text <- gtext(text = "<fist lines>",
                                            container = fiTbl, height = 10))
        size(fiW$text) <- c(600, 100)
        ##===================================================================
        fiTbl[5, 1] <- ( fiW$skipLab <- glabel("Skip lines", container = fiTbl))
        fiTbl[5, 2:2] <- ( fiW$skip <- gspinbutton(from = 0, to = 10, by = 1,
                                                   value = 0, container = fiTbl))
        ##===================================================================
        fiTbl[6, 1] <- ( fiW$headerLab <- glabel("Header?", container =fiTbl) )
        fiTbl[6, 2:6] <- ( fiW$header <- gcheckbox("header ON",
                                                   checked = TRUE,
                                                   container = fiTbl) )
        ##===================================================================
        fiTbl[7, 1] <- ( fiW$decLab <- glabel("decimal sep", container = fiTbl))
        fiTbl[7, 2:6] <- ( fiW$dec <- gradio(items = names(re$RLprov$.decs),
                                             selected = 1,
                                             container = fiTbl,
                                             horizontal = TRUE) )
        ##===================================================================
        fiTbl[8, 1] <- ( fiW$sepLab <- glabel("field sep", container = fiTbl))
        fiTbl[8, 2:6] <- ( fiW$sep <- gradio(items = names(re$RLprov$.seps),
                                             selected = 1,
                                             container = fiTbl,
                                             horizontal = TRUE) )
        ##===================================================================
        ( colTbl <- glayout(container = wFile,
                            label = "Csv file import",
                        cexpand = TRUE, fill = "both") )
        
        colTbl[1, 2] <- (fiW$colLab <- glabel("column", container = colTbl))
        colTbl[1, 3] <- (fiW$formatLab <- glabel("format", container = colTbl))
        colTbl[1, 4] <- (fiW$nameLab <- glabel("name", container = colTbl))
        
        
        colTbl[2, 1] <- (fiW$dateLab <- glabel("date", container = colTbl))
        
        colTbl[2, 2] <- (fiW$dateCol <- gcombobox(items = "<none>                ",
                                                  selected = 1,
                                                  horizontal = TRUE,
                                                  container = colTbl))
        
        colTbl[2, 3] <- (fiW$dateFormat <- glabel("<none>    ", container = colTbl))
        colTbl[2, 4] <- (fiW$dateName <- glabel("date", container = colTbl))
        
        colTbl[3, 1] <- (fiW$dateLab <- glabel("variable", container = colTbl))
        colTbl[3, 2] <- (fiW$varCol <- gcombobox(items = "<none>                ",
                                                 selected = 1,
                                                 horizontal = TRUE,
                                                 container = colTbl))
        colTbl[3, 3] <- (fiW$varFormat <- glabel(" ", container = colTbl))
        colTbl[3, 4] <- (fiW$varName <- gedit("", width = 16, container = colTbl))
        
        ## gseparator(horizontal = TRUE, container = wFile, expand=TRUE)
        ##====================================================================
        (gg <- ggroup(horizontal = TRUE, expand = FALSE,
                      container = wFile))
        addSpring(gg)
        fiW[["Cancel"]] <- gbutton("Cancel", container = gg,
                                   handler = function(h, ...) { dispose(wFile) } )
        
        
        ## item "window" needs to exist prior to the definition of the handler
        fiW[["window"]] <- wFile
        
        fiW[["OK"]] <- gbutton("OK", container = gg,
                               handler = validNewCsvProject,
                               action = fiW)
        ##====================================================================
        fiW$sb <- gstatusbar(paste("Choose a \"fit\" name and edit the eff.",
                                   "duration if needed."),
                             container = wFile)
        ##====================================================================
        ## initialize a flag telling if data can be read or not
        
        enabled(fiW$OK) <- FALSE
        
        ##====================================================================
        ## handlers
        ## A loop would be nice if it worked, even in tcl-tk!!!
        ##====================================================================
        
        addhandlerchanged(obj = fiW$file, handler = tryReadCsv,
                          action = c(fiW, source = "file"))
        addhandlerchanged(obj = fiW$dec, handler = tryReadCsv,
                          action = c(fiW, source = "dec"))
        addhandlerchanged(obj = fiW$sep, handler = tryReadCsv,
                          action = c(fiW, source = "sep"))
        addhandlerchanged(obj = fiW$skip, handler = tryReadCsv,
                          action = c(fiW, source = "skip"))
        addhandlerchanged(obj = fiW$header, handler = tryReadCsv,
                          action = c(fiW, source = "header"))
        addhandlerchanged(obj = fiW$dateCol, handler = tryReadCsv,
                          action = c(fiW, source = "dateCol"))
        addhandlerchanged(obj = fiW$varCol, handler = tryReadCsv,
                          action = c(fiW, source = "varCol"))
        
        addHandlerMouseMotion(obj = fiTbl,
                              handler = function(h, ...)
                                  svalue(fiW$sb) <- "")
        
        ## helps <- c(skip = "number of the lines to skip among the first ones",
        ##           header = "is the first non skipped line a header indicating colnames?")
        
        for (obj in c("skipLab")) {
            hfun <-  function(h, ...) {
                svalue(fiW$sb) <- "number of the lines to skip among the first ones"
            }
            addHandlerMouseMotion(obj = fiW[[obj]], handler = hfun)
        }
        
        for (obj in c("headerLab")) {
            hfun <-  function(h, ...) {
                svalue(fiW$sb) <-
                    "is the first non skipped line a header indicating colnames?"
            }
            addHandlerMouseMotion(obj = fiW[[obj]], handler = hfun)
        }
        for (obj in c("sepLab")) {
            hfun <-  function(h, ...) {
                svalue(fiW$sb) <- "fields (columns) separator in the file"
            }
            addHandlerMouseMotion(obj = fiW[[obj]], handler = hfun)
        }
        
        if (FALSE) {
            for (obj in c("dateColLab")) {
                addHandlerMouseMotion(obj = fiW[[obj]],
                                      handler = function(h, ...)
                                          svalue(fiW$sb) <-
                                              "'date' column: optional but needed for plots")
            }
            
            for (obj in c("dateFormatLab")) {
                addHandlerMouseMotion(obj = fiW[[obj]],
                                      handler = function(h, ...)
                                          svalue(fiW$sb) <- "date format is guessed")
            }
            
            for (obj in c("varCol")) {
                addHandlerMouseMotion(obj = fiW[[obj]],
                                      handler = function(h, ...)
                                          svalue(fiW$sb) <- "'variable' column: REQUIRED")
            }
        }
        
        visible(wFile) <- TRUE
        
    }
    
    
    ##================================================================
    ## Explain how to use an object
    ##  
    ##================================================================           
    
    explain <- function(h, ...) {
        
        svalue(datW$sb) <- h$action
        
    }

    ##================================================================
    ## validate a numerical list
    ##
    ## Here the hander 'h' has
    ##
    ## - 'obj'      the widget that open the file dialog,
    ## - 'action'   the action of the widget (copied to be passed
    ##              to the 'actReadCsvFile' function
    ## 
    ##================================================================           
    
    validNumList <- function(h, ...) {
        
        text <- svalue(h$obj)
        text <- gsub(" ", "", text)
        L <- as.numeric(unlist(strsplit(text, split = ";")))
        
        if (length(L)>0) {
            if (any(is.na(L))) {
                ##alarm()
                svalue(datW$sb) <- "Invalid data. Give numerical values separated by \";\""
            } else {
                svalue(datW$sb) <- sprintf("Numeric list with %d items", length(L))
            }
        } else svalue(datW$sb) <- "Empty list."  
        
    }
    

    ##=========================================================================
    ## delete a project
    ##========================================================================= 
    
    delProject <- function(h, ...) {
        
        if (DEBUG) cat("entering 'delProject'\n")
        
        projName <- svalue(datW$project)
        
        del <- gconfirm(message = sprintf("delete project '%s'?", projName),
                        title = "Delete project")
        
        if (del) {
            i <- svalue(datW$project, index = TRUE)
            ## set another value 
            L <- length(datW$project)
            cat("L =", L, "i =", i, "\n")
            
            if (L > 1) {
                svalue(datW$project) <- datW$project[-i][1]
                
                ## MODIFS 2011-10-14 COMMENT NEXT LINE as in delFit
                ## updateData(h = list(obj = datW$project, action = datW))
                
                ## cat("x", datW$project[1:L], "\n")
                datW$project[] <- datW$project[-i]
                ##cat("y", datW$project[], "\n")      
            } else {
                
                if (TRUE) {
                    ## MODIFS 2011-10-14
                    for ( name in names(re$RLprov$.defaultProj) ) {
                        svalue(datW[[name]]) <- re$RLprov$.defaultProj[[name]]
                        enabled(datW[[name]]) <- FALSE
                    }
                    enabled(datW[["OK"]]) <- FALSE
                    
                    svalue(datW$project, index = TRUE) <- 0
                    datW$project[] <- character(0) 
                } else {     
                    cat("Last item can not be removed")
                    return()
                    cat("last item!!!\n")
                    datW$project[ ] <- "<none>"
                    svalue(datW$project) <- "<none>"
                    svalue(datW$file) <- ""
                    svalue(datW$main) <- ""
                }
            }
            
            ##
            re$dataList[[projName]] <- NULL
            
        }
        
        if (DEBUG) cat("exiting 'delProject'\n")
        
    }
    
    ##======================================================================
    ## Fitting
    ##
    ##======================================================================
    
    fitProj <- function(h, ... ) {
        
        ## print(resW$fit)
        if (DEBUG) cat("entering 'fitProj'\n")
        
        fitName <- svalue(fitW[["fit"]])
        
        if ( !length(fitName) || is.na(fitName) ) {
            alarm()
            svalue(datW$sb) <- "No project selected"
        }
        
        nm <- svalue(fitW$project)
        dL <- re$dataList[[nm]]
        if(DEBUG) cat(sprintf("Fitting with %s\n", nm))
        
        fL <- re$fitsList[[fitName]]
        
        if(DEBUG) print(head(re$dataList[[nm]]$data))
        
        ##conf.pct <- as.numeric(strsplit(svalue(fitW$confLev), ";")[[1]])
        pct.conf <- as.numeric(strsplit(fL$confLev, split = ";")[[1]])
        pred.periods <-  as.numeric(strsplit(fL$returnPeriods, split = ";")[[1]])
        
        if (DEBUG) cat("useMax1 = ", fL$useMax1, "\n")
        if (DEBUG) cat("useMax2 = ", fL$useMax2, "\n")
        
        MAX.data <- NULL
        MAX.effDuration <- NULL
        
        if ( fL$useMax1 == "yes" ) {
            rMax1 <-  as.numeric(unlist(strsplit(dL$rMax1, split = ";")))
            if (DEBUG) cat("rMax1 = ", rMax1, "rmax1Dur = ", dL$rMax1Dur,"\n")
            MAX.data <- c(MAX.data, list(block1 = rMax1))
            MAX.effDuration <- c(MAX.effDuration, as.numeric(dL$rMax1Dur))
        }
        
        if ( fL$useMax2 == "yes" ) {
            rMax2 <-  as.numeric(unlist(strsplit(dL$rMax2, split = ";")))
            if (DEBUG) cat("rMax2 = ", rMax2, "rMax2Dur = ", dL$rMax2Dur,"\n")
            MAX.data <- c(MAX.data, list(block2 = rMax2))
            MAX.effDuration <- c(MAX.effDuration, as.numeric(dL$rMax2Dur))
        }
        
        if (DEBUG) {
            cat("pct.conf = ", pct.conf, "\n")
            cat("pred.periods = ", pred.periods, "\n")
            cat("MAX.data\n")
            print(MAX.data)
            cat("MAX.effDuration\n")
            print(MAX.effDuration)
            cat("Calling 'Renouv'\n")
        }
        
        res <- try(Renouv(x = dL$data$x,
                          effDuration = as.numeric(fL$effDuration),
                          distname.y = as.character(fL$distribution),
                          MAX.data = MAX.data,
                          MAX.effDuration = MAX.effDuration,
                          pct.conf = pct.conf,
                          pred.period = pred.periods,
                          threshold = as.numeric(fL$threshold),
                          main = as.character(fL$main),
                          trace = 0))
        
        ## Add graphical options to the results
        res$rlOptions <- list(defMain = as.character(fL$main),
                              main = as.character(fL$main),
                              xLimSet = "auto",
                              xLimType = "Return period",
                              xLimMin = NA,
                              xLimMax = NA,
                              xLabSet = "auto",
                              defxLab = "years",
                              xLab = "years",
                              yLimSet = "auto",
                              yLimMin = NA,
                              yLimMax = NA,
                              yLabSet = "auto",
                              defyLab = dL$info$varName,
                              yLab = dL$info$varName)
        
        
        if (is(res, "try-error")) {
            
            cat("FIT ERROR\n")
            alarm() 
            svalue(datW$sb) <- geterrmessage()
            
        } else {
            
            if ( fitName %in% names(re$resList) ) {
                
                nms <- names(resW$nB)
                num <- match(x = fitName, table = nms)
                
                if (DEBUG) {
                    cat("fitName = ", fitName, "\n",
                        "nms = ", nms, "\n",
                        "num = ", num,  "\n")
                    cat(sprintf("items in notebook. Item %d will be deleted\n", num))
                    print(nms)
                }
                
                re$resList[[fitName]] <- res
                svalue(resW$nB) <- num        
                
                for (i in 1:length(re$.nB[[fitName]]$children) ) {
                    
                    if (DEBUG) {
                        cat("trying to delete\nwithin ")
                        print(re$.nB[[fitName]]$root)
                        cat("widget\n")
                        print(re$.nB[[fitName]]$children[[i]])
                    }
                    
                    delete(obj = re$.nB[[fitName]]$root,
                           widget = re$.nB[[fitName]]$children[[i]])
                    
                }
                
                if (DEBUG) print(names(resW$nB))
                
                ## dispose(resW[[fitName]]$group)     
                svalue(datW$sb) <- "Existing fit replaced." 
                
            } else {
                
                svalue(datW$sb) <- "New fit added."
                re$resList[[fitName]] <- res
                
                ## copy
                gg1  <- ggroup(label = fitName,
                               horizontal = FALSE,
                               container = resW$nB,
                               expand = TRUE,
                               fill = "both")
                
                re$.nB[[fitName]]$root <- gg1
                
            } 
            
            if (DEBUG) cat("creating in notebook\n")
            
            if (DEBUG) {
                print(names(resW))
                print(resW[[fitName]]$root)
            }
            
            gg2 <- gpanedgroup(label = fitName,
                               horizontal = TRUE,
                               container = re$.nB[[fitName]]$root,
                               expand = TRUE,
                               fill = "x")
            
            text <- format(summary(res))
            
            ## Summary part
            gSum <- gtext(text, container =gg2)
            size(gSum) <- c(600, 280)
            
            gCoef <- gtable(items = makeSummary(results = res),
                            multiple = TRUE,
                            expand = TRUE,
                            fill = "both",
                            container = gg2)
            
            gPred <- gtable(items = roundPred(res$pred),
                            multiple = TRUE,
                            container = re$.nB[[fitName]]$root,
                            expand = TRUE,
                            fill = "both")
            
            re$.nB[[fitName]]$children <- list()
            re$.nB[[fitName]]$children[[1]] <- gg2
            re$.nB[[fitName]]$children[[2]] <- gPred
            
            svalue(datW$sb) <- "New fit added."
            
            ## add3rdMousePopupmenu(gSum, menulist = mL, action = c("summary", fitName))
            
            add3rdMousePopupmenu(gCoef, menulist = menuCopyCoef)
            add3rdMousePopupmenu(gPred, menulist = menuCopyPred)
            
            if (DEBUG) {
                cat("names(re$.nB)\n")
                print(names(re$.nB))
            }
            
        }
        
        if (DEBUG) cat("exiting 'fitProj'\n")
        
    }
    
    
    ##======================================================================
    ## Fitting
    ##
    ##======================================================================
    
    newfitFromProject <- function(h, ...) {
        
        projName <- svalue(h$obj)
        if (DEBUG) cat("project name = ", projName, "\n")
        
        wL <- h$action
        if (DEBUG) cat("widgetList names\n", names(wL), "\n")
        
        fitName <- guessNames(cand = projName, existing = names(re$fitsList))
        svalue(wL$fit) <- fitName
        
        item <- re$dataList[[projName]]
        
        if (DEBUG) cat("names(item) = ", names(item), "\n")
        
        ## XXX
        ## cat("changing thresh. sel.\n")
        ## wL$threshold[ ] <- seq(from = item$thresholds$from,
        ##                       to = item$thresholds$to,
        ##                       by = item$thresholds$by)
        
        svalue(wL$effDuration) <- item$effDuration
        ## svalue(wL$threshold) <- item$thresholds$from
        ## cat("setting 'Main' to ", item$info$shortLab, "\n")
        ## svalue(wL$main) <- item$info$shortLab
        
    }

    ##=====================================================================
    ## Change the definition fit
    ##
    ## using the name of the fit given in 'h$action' assumed to be a
    ## valid name in 'dataList', the widgets are refreshed
    ##=====================================================================
    
    updateFit <- function(h, ...) {
        
        if (DEBUG) cat("[ Entering 'updateFit'\n")
        
        if (DEBUG) {
            cat("values for object 'h'\n")
            print(h$obj[])
        }
        
        if (length(h$obj[])) {
            fitName <- svalue(h$obj)
        } else {
            fitName <- NA
        }
        
        if (DEBUG) cat("fit:     ", fitName, "\n")
        
        if ( length(fitName) && !is.na(fitName) ) {
            
            item <- re$fitsList[[fitName]]
            project <- item$project
            
            if (DEBUG) cat("project: ", project, "\n")
            
            if (DEBUG) cat("   modify widget values: effDuration, ... \n")
            
            fitW$threshold[ ] <- seq(from = re$dataList[[project]]$thresholds$from,
                                     to = re$dataList[[project]]$thresholds$to,
                                     by = re$dataList[[project]]$thresholds$by)
            
            ## for (name in c("project", "effDuration", "useMax1", "useMax2", "threshold", "distribution", "main",
            ##                "confLev", "returnPeriods")) {
            ##   if (DEBUG) cat(sprintf("fitW[[%s]] = %s\n", name, item[[name]]))
            ##   svalue(fitW[[name]]) <- item[[name]]
            ## }
            svalue(fitW[["project"]]) <- item[["project"]]
            svalue(fitW[["effDuration"]]) <- item[["effDuration"]]
            
            svalue(fitW[["useMax1"]]) <- item[["useMax1"]]
            svalue(fitW[["useMax2"]]) <- item[["useMax2"]]
            
            svalue(fitW[["rMax1Lab"]]) <- item[["rMax1Lab"]]
            svalue(fitW[["rMax2Lab"]]) <- item[["rMax2Lab"]]
            
            svalue(fitW[["threshold"]]) <- item[["threshold"]]
            svalue(fitW[["distribution"]]) <- item[["distribution"]]
            svalue(fitW[["confLev"]]) <- item[["confLev"]]
            svalue(fitW[["returnPeriods"]]) <- item[["returnPeriods"]]
            
            enabled(fitW$useMax1) <- ( item[["rMax1Lab"]] != "no value" )
            enabled(fitW$useMax2) <- ( item[["rMax2Lab"]] != "no value" )
            
            enabled(fitW$rMax1Lab) <- ( item[["rMax1Lab"]] != "no value" )
            enabled(fitW$rMax2Lab) <- ( item[["rMax2Lab"]] != "no value" )
            
            for (name in c("effDuration", "threshold", "distribution",
                           "confLev", "returnPeriods")) {
                enabled(fitW[[name]]) <- TRUE
            }
            
            ## svalue(fitW$rMax1Lab) <- dataList[[project]]$rMax1Lab
            ## svalue(fitW$rMax2Lab) <- dataList[[project]]$rMax2Lab   
            
            enabled(fitW$OK) <- TRUE   
            
            if (DEBUG) cat("   ... done\n")
            
        } else {
            
            if (DEBUG) cat("reset\n")
            
            ## resset values to default
            svalue(fitW$project) <- ""
            
            ## quite arbitrary...
            fitW$threshold[] <- seq(from = 0, to = 1, by = 0.1)
            svalue(fitW$threshold) <- 0
            
            svalue(fitW$effDuration) <- ""
            svalue(fitW$useMax1) <- "no"
            svalue(fitW$useMax2) <- "no"
            svalue(fitW$rMax1Lab) <- ""
            svalue(fitW$rMax2Lab) <- ""
            
            enabled(fitW$OK) <- FALSE
            
            if (DEBUG) cat("done\n")
            
            
        }
        
        if (DEBUG) cat("] Exiting 'updateFit'\n")
        
    }
    
    ##===============================================================
    ## Dialog to create and validate a new Fit from a Project
    ##===============================================================
    
    validNewFit <- function(h, ...) {
        
        if (DEBUG) cat("[ Entering 'validNewFit'\n")
        
        ## widgets List
        wL <- h$action
        fitName <- svalue(wL$fit)
        if (DEBUG) cat("fitName = ", fitName, "\n")
        
        ## check duration 
        effDuration <- as.numeric(svalue(wL$effDuration))
        
        if ( (length(effDuration) == 0) || is.na(effDuration) || (effDuration < 0) ) {
            alarm( )
            svalue(wL$sb) <- "Bad effective duration. Must be a positive numeric."
            return()
        }
        
        if ( fitName %in% names(re$fitsList) ) {
            alarm()
            svalue(wL$sb) <- "Fit name already exists!!! Please change it."
            return()
        } else {
            
            if (DEBUG) cat("new name\n")
            
            project <- svalue(wL$project)
            
            ## intialise threshold
            threshold <- re$dataList[[project]]$thresholds$from
            
            ## History 1 and 2
            if ( length(re$dataList[[project]]$rMax1) && nchar(re$dataList[[project]]$rMax1) ) {
                useMax1 <- "yes"
            } else useMax1 <- "no"
            
            if ( length(re$dataList[[project]]$rMax2) && nchar(re$dataList[[project]]$rMax2) ) {
                useMax2 <- "yes"
            } else useMax2 <- "no"
            
            if (DEBUG) {
                cat("useMax1", useMax1, "\n")
                cat("rMax1", re$dataList[[project]]$rMax1, "\n")
                cat("useMax2", useMax2, "\n")
                cat("rMax2", re$dataList[[project]]$rMax2, "\n")
            }
            
            ## enabled(fitW$useMax1) <- ( useMax1 == "yes" )
            ## enabled(fitW$useMax2) <- ( useMax2 == "yes" )
            
            ## svalue(fitW$rMax1) <- abbrNumList(dataList[[project]]$rMax1)
            ## svalue(fitW$rMax2) <- abbrNumList(dataList[[project]]$rMax2)
            
            rMax1Lab <- abbrNumList(re$dataList[[project]]$rMax1)
            rMax2Lab <- abbrNumList(re$dataList[[project]]$rMax2)
            
            ## enabled(fitW$rMax1) <- ( useMax1 == "yes" )
            ## enabled(fitW$rMax2) <- ( useMax2 == "yes" )
            
            newItem <- list(source = "csv file",
                            main = svalue(wL$fit),
                            project = project,
                            threshold = threshold,
                            distribution = "exponential",
                            useMax1 = useMax1,
                            rMax1 = re$dataList[[project]]$rMax1,
                            rMax1Lab = rMax1Lab,
                            useMax2 = useMax2,
                            rMax2 = re$dataList[[project]]$rMax2,
                            rMax2Lab = rMax2Lab,
                            effDuration = svalue(wL$effDuration),
                            confLev = "95;70",
                            returnPeriods = "10;20;50;100;200;500;1000")
            
            svalue(wL$sb) <- "New Fit added"
            if (DEBUG) cat("New fit added\n")
            
            L <- length(re$fitsList)
            
            re$fitsList[[fitName]] <- newItem
            
            ## 2011-12-13 'Lock' the project and disable its gentry
            ## if it is currently selected in the (hidden) tab 'data'
            
            if (DEBUG) cat("Locking the project\n")
            
            re$dataList[[project]]$isLocked <- TRUE
            
            if (svalue(datW$project) == project) {
                updateData(h = list(obj = datW$project, action = datW))
            }
            
            ## end 2011-12-13
            
            ## modify combobox values
            if (L > 0) {
                fitW$fit[1:(L+1)] <- c(fitW$fit[1:L], fitName)
            } else {
                fitW$fit[] <- fitName
            }
            
            if (DEBUG) {
                cat("check\n")
                cat("   L =               ", L, "\n")
                cat("   fitW              ", fitW$fit[ ], "\n")
                cat("   names(re$fitsList)", names(re$fitsList), "\n")
                cat("   fitName           ", fitName, "\n")
                cat("   fitW$fit[]        ", fitW$fit[], "\n")
            }
            
            ## with RGtk2, this works, but not in tcl-tk!
            ## This MUST involve the execution of 'updateFit'
            
            if ( FALSE && options()$guiToolkit == "RGtk2" ) {
                if (DEBUG) cat("svalue(fitW$fit) BEFORE", svalue(fitW$fit), "\n")
                svalue(fitW$fit) <- as.character(fitName)
                if (DEBUG) cat("svalue(fitW$fit) AFTER", svalue(fitW$fit), "\n")
            } 
            
            ## An explicit update was tried in tcl-tk but unsuccessfully
            ##updateFit(h = list(obj = wL$fit, action = fitW))
            
            
        }
        
        ## Later
        dispose(wL[["window"]]) 
        
        if (DEBUG) cat("] Exiting 'validNewFit'\n")
        
    }
  
    ##======================================================================
    ## new fit
    ##
    ##======================================================================

    newFit <- function(h, ...) {
        
        ## fd is for "Fit Dialog"
        fd <- gwindow("New fit", visible = FALSE)
        fdW <- list()
        
        fdTbl <- glayout(container = fd, label = "data")
        
        fdTbl[1, 1] <- glabel("Project", container = fdTbl)
        fdTbl[1, 2] <- ( fdW$project <- gcombobox(items = names(re$dataList),
                                                  selected = 0, container = fdTbl) )
        
        fdTbl[2, 1] <- glabel("Fit", container = fdTbl)
        fdTbl[2, 2] <- ( fdW$fit <- gedit(container = fdTbl) )
        
        fdTbl[3, 1] <- glabel("Eff. duration (yrs)", container =fdTbl)
        fdTbl[3, 2] <- (fdW$effDuration <- gedit(" ", container =fdTbl))
        
        ## The statusbar "sb" is added here because of the evaluation in handler
        fdW$sb <- gstatusbar("Choose a \"fit\" name and edit the eff. duration if needed.",
                             container = fd)
        
        fdW$window <- fd
        
        fdTbl[4, 3] <- (gg <- ggroup(container = fdTbl))
        
        ( fdW$cancel <- gbutton("Cancel", container = gg,
                                handler = function(h, ...) { dispose(fd) } ) )
        
        ( fdW$oK <- gbutton("OK", container = gg,
                            handler = validNewFit,
                            action = fdW) )
        
        addhandlerclicked(obj = fdW$project,
                          handler = newfitFromProject,
                          action = fdW)
        
        visible(fd) <- TRUE
        
    }
    
    ##===========================================================================
    ## delete fit
    ##
    ## XXX Caution: is the deleted fit is the current one, it is necessary to
    ## clean the values of the widgets.
    ##===========================================================================
    
    delFit <- function(h, ...) {
        
        if (DEBUG) cat("[ Entering 'delFit'\n")
        
        fitName <- svalue(fitW$fit)
        
        if(DEBUG) cat("   fit to delete = ", fitName, "\n")
        
        if ( length(fitName) == 0 || is.na(fitName) ) {
            svalue(datW$sb) <- "no active fit to delete!"
            return()
        }
        
        del <- gconfirm(message = sprintf("delete fit '%s'?", fitName),
                        title = "Delete fit")
        
        if (del) {
            
            i <- svalue(fitW$fit, index = TRUE)
            ## set another value 
            L <- length(fitW$fit)
            if (DEBUG) cat("   L =", L, "i =", i, "\n")
            
            if (L > 1) {
                svalue(fitW$fit) <- fitW$fit[-i][1]
                ## updateFit(h = list(obj = fitW$project, action = fitW))
                fitW$fit[] <- fitW$fit[-i]
                
            } else {
                
                if (DEBUG) cat("   last item!!!\n")
                
                ## Changed 2011-10-14 Cleaning fits
                if (FALSE) {
                    svalue(fitW$threshold) <- 0
                    svalue(fitW$effDuration) <- ""
                    svalue(fitW$distribution, index = TRUE) <- 1
                    
                    svalue(fitW$useMax1) <- "no"
                    svalue(fitW$useMax2) <- "no"
                    
                    enabled(fitW$fit) <- FALSE
                    enabled(fitW$useMax1) <- FALSE
                    enabled(fitW$useMax2) <- FALSE
                    
                    
                } else {
                    for ( name in names(re$RLprov$.defaultFit) ) {
                        svalue(fitW[[name]]) <- re$RLprov$.defaultFit[[name]]
                        enabled(fitW[[name]]) <- FALSE
                    }
                }
                
                svalue(fitW$fit, index = TRUE) <- 0
                fitW$fit[] <- character(0)
                svalue(fitW$project) <- ""
                
            }
            
            re$fitsList[[fitName]] <- NULL
            if (DEBUG) cat("   OK\n")
            ## updateFit(h = list(obj = fitW$fit, action = fitW))
        }
        
        if (DEBUG) cat("] Exiting 'delFit'\n")
        
    }
    
    ##=======================================================================
    ## Modify 'fitsList' from the widgets svalue 
    ##=======================================================================
    
    bindFits <- function(h, ...) {
        
        if (DEBUG) cat("[ Entering 'bindFits'\n")
        
        ## MODIFS DU 2011-10-14 NA
        fit <- svalue(fitW$fit)
        if (DEBUG) cat("fit = ", fit, "name = ", name, "\n")
        
        if ( (length(fit) > 0) && !is.na(fit) ) {
            name <- h$action
            
            re$fitsList[[fit]][[name]] <- svalue(fitW[[name]])
        }
        
        if (DEBUG) cat("] Exiting 'bindFits'\n")
        
    }
    
    ##=======================================================================
    ## Modify 'dataList' from the widgets svalue 
    ##=======================================================================
    
    bindAllFits <- function(h, ...) {
        
        if (DEBUG) cat("[ Entering 'bindAllFits'\n")
        
        project <- svalue(fitW$project)
        name <- h$action
        
        re$dataList[[project]][["main"]] <- svalue(fitW[["main"]])
        re$dataList[[project]][["effDuration"]] <- svalue(fitW[["effDuration"]])
        re$dataList[[project]][["rMax1"]] <- svalue(fitW[["rMax1"]])
        re$dataList[[project]][["rMax1Dur"]] <- svalue(fitW[["rMax1Dur"]])
        re$dataList[[project]][["rMax2"]] <- svalue(fitW[["rMax2"]])
        re$dataList[[project]][["rMax2Dur"]] <- svalue(fitW[["rMax2Dur"]])
        
        if (DEBUG) cat("] Exiting 'bindAllFits'\n")
        
    }
    
    ##**********************************************************************
    ## Results part
    ##
    ##**********************************************************************
    
    ##======================================================================
    ## Manage results
    ##
    ##======================================================================
    
    updateRes <- function(h, ... ) {
        
        nm <- svalue(resW[["Fit"]])
        cat("nm = ", nm, "\n")
        
        ## if (length(nm) && !is.na(nm)) { }
        resW$pred[ ] <- re$resList[[nm]]$pred
        
    }

    ##======================================================================
    ## Export results to an 'odt' file
    ##
    ## CODE DELETED see earlier versions!!!
    ##
    ##======================================================================
    
    
    ##======================================================================
    ## Export results to a directory
    ##
    ##======================================================================
    
    actHTMLRes <- function(h, ...) {
        
        HTMLDir <- svalue(h$action$dir)
        
        if (file.exists(HTMLDir)) {
            
            replace <- gconfirm(message = sprintf("replace existing '%s'?", HTMLDir),
                                title = "Replace directory")
            
            if (DEBUG) cat("replace = ", replace, "\n")
            if (!replace) {
                svalue(datW$sb) <- "export canceled"
                return()
            }
            
        } else {
            
            test <- try(dir.create(HTMLDir))
            
            if (is(test, "try-error")) {
                alarm() 
                svalue(datW$sb) <- geterrmessage()
                return()
            } 
            
        }
        
        ## find res name and store HTMLDir to the resList item
        
        L <- length(resW$nB)   
        num <- svalue(resW$nB)
        
        if ( (L == 0) || !length(num) || is.na(num) ) {
            svalue(datW$sb) <- "no results to export"
            return()  
        }
        
        resName <- names(resW$nB)[num]
        re$res <- re$resList[[resName]]
        re$resList[[resName]]$HTMLDir <- HTMLDir
        
        ## copy css steel sheet(s) to HTMLDir
        ## XXX use "file remove" first to clarify the operation
        
        for (media in c("screen", "print")) {
            
            test <- try(file.copy(from = file.path(system.file("reporting", package = "RenextGUI"),
                                      sprintf("reports-%s.css", media)),
                                  to = HTMLDir,
                                  overwrite = TRUE))
            
            if (is(test, "try-error")) {
                alarm() 
                svalue(datW$sb) <- geterrmessage()
                return()
            }
        }
        
        ## write to "index.html" within  HTMLDir
        myFile <- file.path(HTMLDir, "index.html")
        
        cat(paste("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"",
                  "\"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n\n"),
            file = myFile, append = FALSE)
        
        cat(paste("<!-- ========================================================\n",
                  "      file generated by the RenextGUI package for R          \n",
                  "                                                             \n",
                  "      https://gforge.irsn.fr/gf/project/renext.              \n",                
                  "                                                             \n",
                  "     ======================================================== -->\n"),
            file = myFile, append = TRUE)
        
        cat("<html>\n<head>\n", file = myFile, append = TRUE)
        cat(sprintf("   <title>%s</title>\n", resName),
            file = myFile, append = TRUE)
        
        for (media in c("screen", "print")) {
            cat(sprintf(paste("   <link rel=\"stylesheet\" type=\"text/css\" media=\"%s\"",
                              "href=\"./reports-%s.css\"></link>\n"),
                        media, media),
                file = myFile, append = TRUE)
        }
        
        cat("</head>\n\n", file = myFile, append = TRUE)
        cat("<body>\n", file = myFile, append = TRUE)
        
        ## display informations about the fit 
        cat("<div class=\"info\">\n", file = myFile, append = TRUE)
        versions <- installed.packages()[c("Renext", "RenextGUI"), "Version"]
        cat(sprintf(paste("File generated <span class=\"res\">%s</span>",
                          "using <span class=\"res\">Renext_%s</span>",
                          "and <span class=\"res\">RenextGUI_%s</span>\n"),
                    format(Sys.time()), versions[1], versions[2]),
            file = myFile, append = TRUE)
        cat("</div>\n", file = myFile, append = TRUE)
        
        cat("<h2>Estimated parameters</h2>\n",
            file = myFile, append = TRUE)
        
        coefTable <- makeSummary(results = re$res)
        
        HTML(x = coefTable, file = myFile,
             align = "left",
             row.names = FALSE,               ## because par. names are in col #1!
             nsmall = c(3, 5, 2, rep(2, 4)),
             append = TRUE)
        
        cat("<h2>Return levels</h2>\n",
            file = myFile, append = TRUE)
        
        M <- as.matrix(format(roundPred(re$res$pred)))
        rownames(M) <- NULL
        HTML(x = M, file = myFile,
             align = "left",
             row.names = FALSE,
             nsmall = c(3, 5, 2, rep(2, 4)),
             append = TRUE)
        
        PNG <- TRUE
        
        ## Note that we write (later) to HTMLdit, but the link must
        ## be relative since the directory must be movable
        
        if (PNG) {
            graphName <- "RLplot.png"
            png(filename = file.path(HTMLDir, graphName),
                width = 500, height = 400)
        } else {
            graphName <- "RLplot.png"
            svg(filename = file.path(HTMLDir, graphName),
                width = 6, height = 7)
        }
        
        plotRes(h)
        
        dev.off()
        
        R2HTML::HTMLInsertGraph(GraphFileName = file.path(".", graphName),
                                Align = "left",
                                file =  myFile,
                        append = TRUE)
        
        cat("\n\n", file = myFile, append = TRUE)
        
        cat("<h2>Summary</h2>\n",
            file = myFile, append = TRUE)
        
        HTML(x = summary(re$res),
             pred = FALSE, coef = FALSE,
             symbolic.cor = TRUE,
             file = myFile,
             class = "summary",
             append = TRUE)
        
        cat("</body>\n", file = myFile, append = TRUE)
        cat("</html>\n", file = myFile, append = TRUE)
        
        svalue(datW$sb) <- "results have been exported"
        
        
        dispose(h$action$window)
        
        browseURL(paste("file://", myFile, sep = ""))
        
    }
    
    ##======================================================================
    ## Export results to an 'html' directory
    ##
    ##======================================================================
  
    HTMLRes <- function(h, ...) {
        
        if (!requireNamespace("R2HTML", quietly = TRUE)) {
            alarm()
            svalue(datW$sb) <- "package 'R2HTML' not accessible. No HTML export possible!"
            return()
        } 
        
        svalue(datW$sb) <- "HTML export of results"
        
        if (DEBUG) { cat("current item", svalue(resW$nB), "\n") }
        
        L <- length(resW$nB)   
        num <- svalue(resW$nB)
        
        if ( (L > 0) && length(num) && !is.na(num) ) {
            
            resName <- names(resW$nB)[num]
            if (DEBUG) cat(sprintf("Exporting results with name = %s...\n", resName))
            
            if (is.null(cand <- re$resList[[num]]$HTMLDir)) {     
                name <- sub("\\.", "_", resName)
                name <- file.path(getwd(), name)
            } else {
                name <- cand
            }
            
            if (DEBUG) cat("Candidate file name = ", name, "\n")
            
            re$res <- re$resList[[num]]
            
            ## fd is for "Fit Dialog"
            html <- gwindow("HTML export", visible = FALSE)
            htmlW <- list(device = "none")
            
            htmlTbl <- glayout(container = html, label = "data")
            
            htmlTbl[1, 1] <- glabel("Directory", container = htmlTbl)     
            htmlTbl[1, 2:4] <-
                ( htmlW$dir <- gfilebrowse("Dir_to_save",
                                           container = htmlTbl,
                                           initialfilename = name,
                                           width = 50,
                                           type = "selectdir", quote = FALSE) )
            
            ## The statusbar "sb" is added here because of the evaluation in handler
            htmlW$sb <- gstatusbar("Choose an export directory.",
                                   container = html)
            
            htmlW$window <- html
            
            htmlTbl[4, 3] <-
                ( htmlW$cancel <- gbutton("Cancel",
                                          container = htmlTbl,
                                          handler = function(h, ...) { dispose(html) } ) )
            
            htmlTbl[4, 4] <-
                ( htmlW$oK <- gbutton("OK", container = htmlTbl,
                                      handler = actHTMLRes,
                                      action = htmlW) )
            
            visible(html) <- TRUE
            
            
        } else {
            alarm()
            svalue(datW$sb) <- "no results to export"
        }
        
    }
    
    ##======================================================================
    ## Actually dump the fit as a R program
    ##
    ##======================================================================
    
    actDumpRes <- function(h, ...) {
        
        fileName <- svalue(h$action$file)
        
        test <- try(dumpFitAndRes(resName = h$action$resName,
                                  fileName = fileName))
        
        if (is(test, "try-error")) {
            alarm() 
            svalue(datW$sb) <- geterrmessage()
            return()
        }
        
        dispose(h$action$window)
        
        ## Next line removed because it creates problems under Windows 
        ## browseURL(paste("file://", fileName, sep = ""))
        
        
    }
    
    ##======================================================================
    ## Dialog before dumping the fit as a R program
    ##
    ##======================================================================
    
    dumpRes <- function(h, ...) {
        
        
        svalue(datW$sb) <- "R 'dump' of the fit/results"
        
        if (DEBUG) { cat("current item", svalue(resW$nB), "\n") }
        
        L <- length(resW$nB)   
        num <- svalue(resW$nB)
        
        if ( (L > 0) && length(num) && !is.na(num) ) {
            
            resName <- names(resW$nB)[num]
            if (DEBUG) cat(sprintf("Dumping fit/ results with name = %s...\n", resName))
            
            if (is.null(cand <- re$resList[[num]]$Rprog)) {     
                name <- sub("\\.", "_", resName)
                name <- file.path(getwd(), name)
            } else {
                name <- cand
            }
            
            if (DEBUG) cat("Candidate file name = ", name, "\n")
            
            ## fd is for "Fit Dialog"
            Rdump <- gwindow("R dump", visible = FALSE)
            rW <- list(resName = resName)
            
            rTbl <- glayout(container = Rdump, label = "data")
            
            rTbl[1, 1] <- glabel("File", container = rTbl)     
            rTbl[1, 2:4] <-
                ( rW$file <- gfilebrowse("File_to_save",
                                         container = rTbl,
                                         initialfilename = name,
                                         filter = list("R files" = list(patterns = c("*.[rR]"))),
                                         width = 50,
                                         type = "save",
                                         quote = FALSE) )
            
            ## The statusbar "sb" is added here because of the evaluation in handler
            rW$sb <- gstatusbar("Choose a file",
                                container = Rdump)
            
            rW$window <- Rdump
            
            rTbl[4, 3] <-
                ( rW$cancel <- gbutton("Cancel",
                                       container = rTbl,
                                       handler = function(h, ...) { dispose(Rdump) } ) )
            
            rTbl[4, 4] <-
                ( rW$oK <- gbutton("OK", container = rTbl,
                                   handler = actDumpRes,
                                   action = rW) )
            
            visible(Rdump) <- TRUE
            
            
        } else {
            alarm()
            svalue(datW$sb) <- "no fit/results to dump"
        }
        
    }
    
    ##======================================================================
    ## Dump a result as a R code using the result name to find fit and
    ## project properties
    ##
    ## if eval is TRUE, the n the code will be evaluated.
    ##
    ##======================================================================

    dumpFitAndRes <- function(resName, fileName = NULL, eval = FALSE) {
        
        fit <- re$fitsList[[resName]]
        projName <- fit$project
        proj <- re$dataList[[projName]]
        res <- re$resList[[resName]]
        
        version <- installed.packages()[c("RenextGUI"), "Version"]
        
        myText <-
            paste("##*****************************************************************************",
                  "## File generated on %s by RenextGUI version %s ",
                  "## ",
                  "## CAUTION: some objects from the global environment are replaced, such as ",
                  "## 'myData', 'myFit',... ",
                  "## ",
                  "## You might change the name of the created objects 'myData', 'myFit'.",
                  "## ",
                  "##*****************************************************************************\n",
                  sep = "\n")
        
        myCode <- sprintf(fmt = myText, format(Sys.time()), version)
        
        myText <- "library(Renext)\n"
        
        myCode <- paste(myCode, myText, sep = "\n")
        
        myText <-
            paste("##=============================================================================",
                  "## Read data from the csv file, using the parameters provided by the GUI       ",
                  "##                                                                             ",
                  "## Setting 'stringAsFactor' to FALSE avoids the implicit conversion of         ",
                  "## strings to factors.                                                         ",
                  "##                                                                             ",
                  "## Note that using 'colClass' arg, we would have read only the wanted columns  ",
                  "## in a classical programming context.                                         ",
                  "##=============================================================================\n",
                  sep = "\n")
        
        myCode <- paste(myCode, myText, sep = "\n")
        
        myText <-
            paste("myData <- ",
                  "    read.table(file = \"%s\",",
                  "               sep = \"%s\",",
                  "               header = %s,",
                  "               skip = %s,",
                  "               dec = \"%s\",",
                  "               stringsAsFactors = FALSE)\n",
                  sep = "\n")
        
        myCode <- paste(myCode,
                        sprintf(fmt = myText,
                                proj$path, proj$csvPar$sep,
                                proj$csvPar$header, proj$csvPar$skip,
                                proj$csvPar$dec),
                        sep = "\n")
        
        ##====================================================================
        ## generate the suitable code, according the existence of a date
        ##====================================================================
        
        if (proj$hasDate) {
            myText <-
                paste("myData <- ",
                      "    data.frame(date = as.POSIXct(strptime(myData[ , %d], format = \"%s\")),",
                      "               x = myData[ , %d])",
                      "",
                      sep = "\n")
            
            myText <- sprintf(myText,
                              as.integer(proj$colDate),
                              proj$dateFormat,
                              as.integer(proj$colNum))
        } else {
            myText <- "myData <- data.frame(x = myData[ , %d])\n"
            myText <- sprintf(myText, as.integer(proj$colNum))
        }
        
        myCode <- paste(myCode, myText, sep = "\n")
        
        myCode <-
            paste(myCode,
                  paste("##=============================================================================",
                        "## fit a Renouv model using the parameters defined through the GUI             ",
                        "## Caution: some rounding may occur for numeric values                         ",
                        "##=============================================================================\n",
                        sep = "\n"),
                  sep = "\n")
        
        predPeriods <-  paste("c(",
                              paste(strsplit(fit$returnPeriods, split = ";")[[1]], collapse = ", "),
                              ")", sep = "")
        
        confLev <-  paste("c(",
                          paste(strsplit(fit$confLev, split = ";")[[1]], collapse = ", "),
                          ")", sep = "")
        
        MAX.data.str <- NULL
        MAX.effDuration.str <- NULL
        
        if ( fit$useMax1 == "yes" ) {
            rMax1.str <-  paste("block1 = c(",
                                paste(unlist(strsplit(proj$rMax1, split = ";")),
                                      collapse = ", "),
                                ")", sep = "")
            ## MAX.data.str <- list(block1 = rMax1)
            MAX.effDuration.str <- c(MAX.effDuration.str, proj$rMax1Dur)
        } else {
            rMax1.str <- character(0)
        }
        
        if ( fit$useMax2 == "yes" ) {
            rMax2.str <-  paste("block2 = c(",
                                paste(unlist(strsplit(proj$rMax2, split = ";")),
                                      collapse = ", "),
                                ")", sep = "")
            ## MAX.data <- c(MAX.data, list(block2 = rMax2))
            MAX.effDuration.str <- c(MAX.effDuration.str, proj$rMax2Dur)
        } else {
            rMax2.str <- character(0)
        }
        
        if ( (fit$useMax1 == "yes") || (fit$useMax2 == "yes") ) {
            
            sep <- ifelse( (fit$useMax1 == "yes") && (fit$useMax2 == "yes"),
                          ", \n     ", "\n")
            
            MAX.data.str <- paste("MAX.data <- list(",
                                  paste(rMax1.str, rMax2.str, sep = sep),
                                  ")", sep  = "\n     ")
            
            MAX.effDuration.str <- paste("MAX.effDuration <- c(",
                                         paste(MAX.effDuration.str, collapse = ", "),
                                         ")", sep  = "")
        } else {
            MAX.data.str <- "MAX.data <- NULL"
            MAX.effDuration.str <- "MAX.effDuration <- NULL"
        }
        
        
        myCode <-
            paste(myCode,
                  MAX.data.str,
                  "",
                  sep = "\n")
        
        myCode <-
            paste(myCode,
                  MAX.effDuration.str,
                  "",
                  sep = "\n")
        
        myCode <-
            paste(myCode, 
                  sprintf(paste("myFit <- ",
                                "    Renouv(x = myData$x,",
                                "           effDuration = %s,",
                                "           distname.y = \"%s\",",
                                "           MAX.data = %s,",
                                "           MAX.effDuration = %s,",
                                "           pct.conf = %s,",
                                "           pred.period = %s,",
                                "           threshold = %s,",
                                "           main = \"%s\",",
                                "           trace = 0)\n",
                                sep = "\n"),
                          fit$effDuration,
                          fit$distribution,
                          "MAX.data",
                          "MAX.effDuration",
                          confLev,
                          predPeriods,
                          fit$threshold,
                          fit$main),
                  sep = "\n")
        
        myText <-
            paste("##=============================================================================",
                  "## Use suitable methods on the result 'myFit', an object with class \"Renouv\" ",
                  "##                                                                             ",
                  "##                 summary, plot, predict, ...                                 ",
                  "##                                                                             ",
                  "## See the documentation of the Renext package                                 ",
                  "##                                                                             ",
                  "##=============================================================================\n",
                  sep = "\n")
        
        myCode <- paste(myCode, myText, sep = "\n")
       
        myCode <-
            paste(myCode,
                  "summary(myFit)",
                  "plot(myFit)",
                  "predict(myFit)",
                  sep = "\n")
        
        ##============================================================================
        ## now, the R code is ready for parsing and evaluation
        ##============================================================================
        
        if (!is.null(fileName)) {
            
            test <- writeLines(text = myCode,
                               con = fileName,
                               sep = "\n", useBytes = FALSE)    
            
        }
        
        ##=============================================================================
        ## 
        ##=============================================================================
        
        if (eval) {
            eval(parse(text = myCode))
        }
        
    }
    
    
    ##======================================================================
    ## copy the table of coefficients for the active results to clipboard 
    ##
    ##======================================================================
    
    coef2Clipboard <- function(h, ...) {
        
        svalue(datW$sb) <- "copy to clipboard"
        
        if (DEBUG) { cat("current item", svalue(resW$nB), "\n") }
        
        L <- length(resW$nB)   
        num <- svalue(resW$nB)
        
        if ( (L > 0) && length(num) && !is.na(num) ) {
            
            resName <- names(resW$nB)[num]
            if (DEBUG) cat(sprintf("Copy to cliboard coef for %s...\n", resName))
            
            res <- try(write.table(x = makeSummary(results = re$resList[[resName]]),
                                   file = "clipboard",
                                   row.names = FALSE,
                                   sep = "\t"))
            if (is(res, "try-error")) {
                alarm() 
                svalue(datW$sb) <- geterrmessage()
            }
            
            if (DEBUG) print(makeSummary(results = re$resList[[resName]]))
            
        }
        
    }
    
    ##======================================================================
    ## copy the table of predictions for the active results to clipboard 
    ##
    ##======================================================================
    
    pred2Clipboard <- function(h, ...) {
        
        svalue(datW$sb) <- "copy to clipboard"
        
        if (DEBUG) { cat("current item", svalue(resW$nB), "\n") }
        
        L <- length(resW$nB)   
        num <- svalue(resW$nB)
        
        if ( (L > 0) && length(num) && !is.na(num) ) {
            
            
            ##=======================
            ## 2011-10-03 
            ## write.table(as.matrix(format(roundPred(resList[[resName]]$pred))))
            ## write.table(x = as.matrix(format(roundPred(resList[[resName]]$pred))),
            ##     file = "clipboard",
            ##     quote = FALSE, sep = "\t")
            
            resName <- names(resW$nB)[num]
            if (DEBUG) cat(sprintf("Copy to clipboard pred for  %s...\n", resName))
            pred <- roundPred(re$resList[[resName]]$pred)
            pred$period <- round(pred$period, digits = 0)
            
            res <- try(write.table(x = pred,
                                   file = "clipboard",
                                   row.names = FALSE,
                                   sep = "\t"))
            
            if (is(res, "try-error")) {
                alarm() 
                svalue(datW$sb) <- geterrmessage()
            }
            if (DEBUG) print(roundPred(re$resList[[resName]]$pred))
            
            
        }
        
    }
    
    ##======================================================================
    ## Plot results
    ##
    ##======================================================================
    
    plotRes <- function(h, ...) {
        
        if (DEBUG) cat("entering 'plotRes'\n")
        svalue(datW$sb) <- "return level plot"
        
        L <- length(resW$nB)   
        num <- svalue(resW$nB)
        
        if ( (L > 0) && length(num) && !is.na(num) ) {
            
            resName <- names(resW$nB)[num]
            if (DEBUG) cat(sprintf("Plotting results with name = %s...\n", resName))
            
            rlOpt <- re$resList[[resName]]$rlOptions
            
            if (rlOpt$xLimSet == "set") {
                if (rlOpt$xLimType == "prob.") {
                    problim <- as.numeric(c(rlOpt$xLimMin, rlOpt$xLimMax))
                    Tlim <- NULL
                } else {
                    Tlim <- as.numeric(c(rlOpt$xLimMin, rlOpt$xLimMax))
                    problim <- NULL
                }
            } else {
                Tlim <- NULL
                problim <- NULL
            }
            
            if (DEBUG) cat("rlOpt$defxLab",  rlOpt$defxLab, "\n")
            
            if (rlOpt$xLabSet == "set") xLab <- rlOpt$xLab
            else xLab <- rlOpt$defxLab
            
            if (rlOpt$yLabSet == "set") yLab <- rlOpt$yLab
            else yLab <- rlOpt$defyLab
            
            ## This could be generalized in the future to
            ## actions before plotting and actions after...
            
            if (h$action$device == "x11") dev.new()
            
            if (rlOpt$yLimSet == "auto") {
                
                ## 'label' must be set for recent versions
                if (packageVersion("Renext") >= "2.0-2") {
                    
                    plot(re$resList[[resName]],
                         Tlim = Tlim, problim = problim,
                         main = rlOpt$main,
                         label = "",
                         xlab = xLab, ylab = yLab)
                    
                } else {
                    
                    plot(re$resList[[resName]],
                         Tlim = Tlim, problim = problim,
                         main = rlOpt$main,
                         xlab = xLab, ylab = yLab)
                    
                }
                
            } else {
                ## cat ("ylim = ", as.numeric(c(rlOpt$yLimMin, rlOpt$yLimMax)), "\n")
                if (packageVersion("Renext") >= "2.0-2") {
                    
                    plot(re$resList[[resName]],
                         Tlim = Tlim, problim = problim,
                         ylim = as.numeric(c(rlOpt$yLimMin, rlOpt$yLimMax)),
                         main = rlOpt$main,
                         label = "",
                         xlab = xLab, ylab = yLab)
                    
                } else {
                    
                    plot(re$resList[[resName]],
                         Tlim = Tlim, problim = problim,
                         ylim = as.numeric(c(rlOpt$yLimMin, rlOpt$yLimMax)),
                         main = rlOpt$main,
                         xlab = xLab, ylab = yLab)
                    
                }
                
            }
            
        } else {
            alarm()
            svalue(datW$sb) <- "no results to plot"
            
        }
        
        if (DEBUG) cat("exiting 'plotRes'\n")
        
    }
    
    
    ##=========================================================================
    ## delete Results
    ##
    ##=========================================================================
    
    delRes <- function(h, ...) {
        
        if (DEBUG) cat("[ Entering 'delRes'\n")
        
        L <- length(resW$nB)    
        num <- svalue(resW$nB)
        
        if ( (L > 0) &&  length(num) && !is.na(num) ) {
            
            resName <- names(resW$nB)[num]
            
            del <- gconfirm(message = sprintf("delete results '%s'?", resName),
                            title = "Delete results")
            
        } else {
            del <- FALSE
            alarm()
            svalue(datW$sb) <- "delete impossible"
            return()
        }
        
        if (del) {
            
            if (DEBUG) cat("note book content\n")
            dispose(resW$nB)
            
            re$resList[[resName]] <- NULL
            svalue(datW$sb) <- "results deleted"
            
            if (DEBUG) cat("OK\n")
        } else {
            svalue(datW$sb) <- ""
        }
        
        if (DEBUG) cat("] Exiting 'delRes'\n")
        
    }
    
    ##=======================================================================
    ## Check and validate options
    ## 
    ##
    ##=======================================================================
    
    validRlOptions <- function(h, ...) {
        
        if (DEBUG) cat("entering 'validRlOptions'\n")
        
        OK <- TRUE
        
        L <- h$action
        newL <- list()
        
        for (nm in c("defMain", "defxLab", "defyLab")) {
            newL[[nm]] <- L[[nm]]
        }
        
        for (obj in c("main",
                      "xLimSet", "xLimType", "xLimMin", "xLimMax",
                      "xLabSet", "xLab",
                      "yLimSet", "yLimMin", "yLimMax",
                      "yLabSet", "yLab")) {
            ## cat(obj, " value = ", svalue(L[[obj]]), "\n")
            newL[[obj]] <- svalue(L[[obj]])
        }
        
        if (newL$xLimSet == "set") {
            
            for ( obj in c("xLimMin", "xLimMax") ) {
                test <- try(x <- as.numeric(newL[[obj]]))
                if ( is(test, "try-error") || is.na(x) ) {
                    alarm(); svalue(L$sb) <- paste("error", obj)
                    OK <- FALSE
                }
            }
            
            if ( (newL$xLimType == "prob.") &&
                ( (as.numeric(newL$xLimMax) > 1.0) || ( (as.numeric(newL$xLimMin) < 0.0) ) )) {
                alarm(); svalue(L$sb) <- paste("invalid probability")
                OK <- FALSE
            } 
            if ( (newL$xLimType == "prob.") && (as.numeric(newL$xLimMin) < 0.0) ) {
                alarm(); svalue(L$sb) <- paste("invalid return period")
                OK <- FALSE
            } 
            if ( as.numeric(newL$xLimMax) <= as.numeric(newL$xLimMin) ) { 
                alarm(); svalue(L$sb) <- paste("'x max' must be greater than 'x min'")
                OK <- FALSE
            }
            
        }
        
        if (newL$yLimSet == "set") {
            
            for ( obj in c("yLimMin", "yLimMax") ) {
                test <- try(x <- as.numeric(newL[[obj]]))
                if (is(test, "try-error") || is.na(x)) {
                    alarm(); svalue(L$sb) <- paste("error", obj)
                    OK <- FALSE
                }
            }
            
            if ( as.numeric(newL$yLimMax) <= as.numeric(newL$yLimMin) ) { 
                alarm(); svalue(L$sb) <- paste("'y max' must be greater than 'y min'")
                OK <- FALSE
            }
            
        }
        
        if (OK) {
            re$resList[[L$name]]$rlOptions <- newL
            dispose(L$window)
        }
        
        if (DEBUG) cat("exiting 'validRlOptions'\n")
        
    }
    

    ##=======================================================================
    ## Dialog for return level options
    ## 
    ##
    ##=======================================================================

    rlOptions <- function(h, ...) { 
        
        
        if (DEBUG) cat("entering 'rlOptions'\n")
        svalue(datW$sb) <- "edit return level options"
        
        if (DEBUG) { cat("current item", svalue(resW$nB), "\n") }
        
        L <- length(resW$nB)   
        num <- svalue(resW$nB)
        
        if ( (L == 0) || (length(num) == 0) || is.na(num) ) {
            svalue(datW$sb) <- "no result seected for options"
            return()
            
        }
        
        rlW <- list( )
        
        rlW$name <- names(resW$nB)[num]
        L <- re$resList[[num]]$rlOptions
        
        for (nm in c("defMain", "defxLab", "defyLab")) {
            rlW[[nm]] <- L[[nm]]
        }
        
        rlOptions <- gwindow("Options for Return Level Plot",
                             height = 500, width = 500, visible = FALSE)
        
        ( rlW$sb <- gstatusbar("edit options for the return level plot",
                               container = rlOptions) )
        
        rlTbl <- glayout(container = rlOptions,
                         label = "Options for Return Level Plot")
        ##---------------------------------------------------------------------------
        rlTbl[1, 1:4] <- gseparator(horizontal = TRUE, container = rlTbl, expand = TRUE)
        ##---------------------------------------------------------------------------
        rlTbl[2, 1] <- glabel("Title", container = rlTbl)
        rlTbl[2, 2] <- ( rlW$main <- gedit(L$main, quote = FALSE, container = rlTbl) )
        ##---------------------------------------------------------------------------
        rlTbl[3, 1:4] <- gseparator(horizontal = TRUE, container = rlTbl, expand = TRUE)
        ##---------------------------------------------------------------------------
        rlTbl[4, 1] <- glabel("x axis limits", container = rlTbl)
        rlTbl[4, 2:3] <- ( rlW$xLimSet  <- gradio(items = c("auto", "set"),
                                                  selected = 1,
                                                  container = rlTbl,
                                                  horizontal = TRUE) )
        ##---------------------------------------------------------------------------
        rlTbl[5, 1] <- ( rlW$xLimTypeLab <-  glabel("limits in...", container =rlTbl) )
        rlTbl[5, 2:3] <- ( rlW$xLimType  <- gradio(items = c("Return period", "prob."),
                                                   selected = 1,
                                                   container = rlTbl,
                                                   horizontal = TRUE) )
        enabled(rlW$xLimTypeLab) <- FALSE
        enabled(rlW$xLimType) <- FALSE
        ##---------------------------------------------------------------------------
        rlTbl[6, 2:5] <- (gx <- ggroup(horizontal = TRUE, expand = FALSE,
                                       container = rlTbl))
        ( rlW$xLimMinLab <- glabel("min", container =gx) )
        ( rlW$xLimMin  <- gedit(L$xLimMin, container =gx) )
        ( rlW$xLimMaxLab <- glabel("max", container =gx)  )
        ( rlW$xLimMax <- gedit(L$xLimMax, container =gx) )
        enabled(rlW$xLimMinLab) <- FALSE
        enabled(rlW$xLimMin) <- FALSE
        enabled(rlW$xLimMaxLab) <- FALSE
        enabled(rlW$xLimMax) <- FALSE
        ##---------------------------------------------------------------------------
        rlTbl[7, 1:4] <- gseparator(horizontal = TRUE, container = rlTbl, expand = TRUE)
        ##---------------------------------------------------------------------------
        rlTbl[8, 1] <- glabel("x label", container = rlTbl)
        rlTbl[8, 2:3] <- ( rlW$xLabSet  <- gradio(items = c("auto", "set"),
                                                  selected = 1,
                                                  container = rlTbl,
                                                  horizontal = TRUE) )
        rlTbl[9, 2] <- ( rlW$xLab  <- gedit(L$xLab, container = rlTbl) )
        enabled(rlW$xLab) <- FALSE
        ##---------------------------------------------------------------------------
        rlTbl[10, 1:4] <- gseparator(horizontal = TRUE, container = rlTbl, expand = TRUE)
        ##---------------------------------------------------------------------------
        rlTbl[11, 1] <- glabel("y axis limits", container =rlTbl)
        rlTbl[11, 2:3] <- ( rlW$yLimSet  <- gradio(items = c("auto", "set"),
                                                   selected = 1,
                                                   container =rlTbl,
                                                   horizontal = TRUE) )
        ##---------------------------------------------------------------------------
        rlTbl[12, 2:5] <- (gy <- ggroup(horizontal = TRUE, expand = FALSE,
                                        container = rlTbl))
        ( rlW$yLimMinLab <- glabel("min", container = gy) )
        ( rlW$yLimMin  <- gedit(L$yLimMin, container = gy) )
        ( rlW$yLimMaxLab <- glabel("max", container = gy) )
        ( rlW$yLimMax <- gedit(L$yLimMax, container = gy) )
        
        enabled(rlW$yLimMinLab) <- FALSE
        enabled(rlW$yLimMin) <- FALSE
        enabled(rlW$yLimMaxLab) <- FALSE
        enabled(rlW$yLimMax) <- FALSE
        ##---------------------------------------------------------------------------
        rlTbl[13, 1:4] <- gseparator(horizontal = TRUE, container =rlTbl, expand = TRUE)
        ##---------------------------------------------------------------------------
        rlTbl[14, 1] <- glabel("y label", container = rlTbl)
        rlTbl[14, 2:3] <- ( rlW$yLabSet  <- gradio(items = c("auto", "set"),
                                                   selected = 1,
                                                   container = rlTbl,
                                                   horizontal = TRUE) )
        rlTbl[15, 2] <- ( rlW$yLab  <- gedit(L$yLab, container = rlTbl) )
        enabled(rlW$yLab) <- FALSE
        ##---------------------------------------------------------------------------
        rlTbl[16, 1:4] <- gseparator(horizontal = TRUE, container =rlTbl, expand = TRUE)
        ##---------------------------------------------------------------------------
        rlTbl[17, 4] <- (gg <- ggroup(horizontal = TRUE, expand = FALSE,
                                      container = rlTbl))
        
        rlW$window <- rlOptions
        
        rlW[["Cancel"]] <- gbutton("Cancel", container =gg,
                                   handler = function(h, ...) { dispose(rlOptions) } )
        
        rlW[["OK"]] <- gbutton("OK", container =gg,
                               handler = validRlOptions,
                               action = rlW)
        
        ##---------------------------------------------------------------------------
        addHandlerMouseMotion(obj = rlTbl,
                              handler = function(h, ...) {svalue(rlW$sb) <- "" } )
        
        addHandlerMouseMotion(obj = rlW[["main"]],
                              handler = function(h, ...) svalue(rlW$sb) <- "give main title")
        
        addHandlerMouseMotion(obj = rlW[["xLimMin"]],
                              handler = function(h, ...) svalue(rlW$sb) <- "give min value")
        
        addHandlerMouseMotion(obj = rlW[["yLimMin"]],
                              handler = function(h, ...) svalue(rlW$sb) <- "give min value")
        
        addHandlerMouseMotion(obj = rlW[["xLimMax"]],
                              handler = function(h, ...) svalue(rlW$sb) <- "give max value")
        
        addHandlerMouseMotion(obj = rlW[["yLimMax"]],
                              handler = function(h, ...) svalue(rlW$sb) <- "give max value")
        
        addHandlerMouseMotion(obj = rlW[["xLab"]],
                              handler = function(h, ...) svalue(rlW$sb) <- "give axis label")
        
        addHandlerMouseMotion(obj = rlW[["yLab"]],
                              handler = function(h, ...) svalue(rlW$sb) <- "give axis label")
        
        ## DOES NOT WORK in tcltk
        ## addHandlerMouseMotion(obj = gg,
        ##                      handler = function(h, ...) svalue(rlW$sb) <- "")
        
        ##============================================================================
        ## handlers
        ##============================================================================
        
        addhandlerchanged(obj = rlW$xLimSet,
                          handler = function(h, ....) {
                              if (svalue(h$obj) == "set") {
                                  enabled(rlW$xLimTypeLab) <- TRUE
                                  enabled(rlW$xLimType) <- TRUE
                                  enabled(rlW$xLimMinLab) <- TRUE
                                  enabled(rlW$xLimMin) <- TRUE
                                  enabled(rlW$xLimMaxLab) <- TRUE
                                  enabled(rlW$xLimMax) <- TRUE
                              } else {
                                  enabled(rlW$xLimTypeLab) <- FALSE
                                  enabled(rlW$xLimType) <- FALSE
                                  enabled(rlW$xLimMinLab) <- FALSE
                                  enabled(rlW$xLimMin) <- FALSE
                                  enabled(rlW$xLimMaxLab) <- FALSE
                                  enabled(rlW$xLimMax) <- FALSE
                              }
                          })
        
        addhandlerchanged(obj = rlW$yLimSet,
                          handler = function(h, ....) {
                              if (svalue(h$obj) == "set") {
                                  enabled(rlW$yLimMinLab) <- TRUE
                                  enabled(rlW$yLimMin) <- TRUE
                                  enabled(rlW$yLimMaxLab) <- TRUE
                                  enabled(rlW$yLimMax) <- TRUE
                              } else {
                                  enabled(rlW$yLimMinLab) <- FALSE
                                  enabled(rlW$yLimMin) <- FALSE
                                  enabled(rlW$yLimMaxLab) <- FALSE
                                  enabled(rlW$yLimMax) <- FALSE
                              }
                          })
        
        addhandlerchanged(obj = rlW$xLabSet,
                          handler = function(h, ....) {
                              if (svalue(h$obj) == "set") {
                                  enabled(rlW$xLab) <- TRUE
                              } else {
                                  enabled(rlW$xLab) <- FALSE
                              }
                          })
        
        addhandlerchanged(obj = rlW$yLabSet,
                          handler = function(h, ....) {
                              if (svalue(h$obj) == "set") {
                                  enabled(rlW$yLab) <- TRUE
                              } else {
                                  enabled(rlW$yLab) <- FALSE
                              }
                          })
        
        ## now set values
        svalue(rlW$xLimSet) <- L$xLimSet
        svalue(rlW$xLimType) <- L$xLimType
        svalue(rlW$xLabSet) <- L$xLabSet
        svalue(rlW$yLimSet) <- L$yLimSet
        svalue(rlW$yLabSet) <- L$yLabSet
        
        visible(rlOptions) <- TRUE
        
    }
    

    ##**************************************************************************
    ## Construction of the widgets of this frame.
    ## 
    ## 
    ##**************************************************************************
    
    
    ##======================================================================
    ## Start the GUI and buit 
    ##
    ##======================================================================
    
    if (DEBUG) cat("Build window and data widgets...\n")
    
    renGUI <- gwindow("Renext GUI", width = 900, visible = FALSE)
    
    ## Does not work correclty at the time
    if (FALSE) {
        tbl <- list(new = list(handler = function(h, ...) { }),
                    delete = list(handler = function(h, ...) { }),
                    plot = list(handler = function(h, ...) { }),
                    quit = list(handler = function(h,...) dispose(renGUI)))
        
        tb <-  gtoolbar(toolbarlist = tbl, container = renGUI)  
    }
    
    noteBook <- gnotebook(container = renGUI, expand = TRUE)
    
    ## maybe add a container here???
    dataTbl <- glayout(container = noteBook, label = "data", expand = TRUE, fill = "both")
    fitTbl <- glayout(container = noteBook, label = "fits")
    resGrp <- ggroup(container = noteBook,  horizontal = FALSE,
                     label = "results", expand = TRUE, fill = "both")
    
    logTbl <- glayout(container = noteBook, label = "log")
    
    ##======================================================================
    ## data part
    ##======================================================================
    
    datW <- list()
    
    ## Buttons group
    dataTbl[1, 1:2] <- ( gd0 <- ggroup(container =dataTbl) )
    gbutton(text = "New", container = gd0, handler = newCsvProject)
    gbutton(text = "Delete", container = gd0, handler = delProject)
    gbutton(text = "Plot", active = FALSE, container = gd0,
            handler = plotData)
    
    ## Project combo
    dataTbl[2:3, 1] <- glabel("Project", container = dataTbl)
    dataTbl[2:3, 2] <- ( datW$project <- gcombobox(items = names(re$dataList),
                                                   selected = 0, container = dataTbl) )
    
    ## File name
    dataTbl[4, 1] <- glabel("OT data csv file", container = dataTbl)
    dataTbl[4, 2] <- (datW$file  <- glabel("                                    ",
                                           container = dataTbl,
                                           ## handler = actReadCsvFile,  ## LATER
                                           action = datW))
    
    ## dataTbl[4, 3] <- gbutton(text = "edit", handler = readCsvFile, container = dataTbl)
    
    ## varName 
    dataTbl[5, 1 ] <- glabel("Var. name", container = dataTbl)
    dataTbl[5, 2 ] <- (datW$varName <- glabel(" ", container = dataTbl))
    
    ## Main 
    dataTbl[6, 1 ] <- glabel("Title", container = dataTbl)
    dataTbl[6, 2 ] <- (datW$main <- gedit(" ", width = 50,
                                          init_msg = "<title>",
                                          container = dataTbl))
    
    dataTbl[7, 1:2] <- gseparator(horizontal = TRUE, container = dataTbl, expand = TRUE)
    ## Effective duration
    dataTbl[8, 1 ] <- glabel("Eff. Dur.", container = dataTbl)
    dataTbl[8, 2 ] <- (datW$effDuration <- gedit(" ", width = 50, container = dataTbl))
    
    dataTbl[9, 1:2] <- gseparator(horizontal = TRUE, container = dataTbl, expand = TRUE)
    
    ## rMax1 data 
    dataTbl[10, 1 ] <- glabel("Hist. Max 1", container = dataTbl)
    dataTbl[10, 2 ] <- (datW$rMax1 <- gedit(" ", width = 50, container = dataTbl))
    
    ## rMax1 duration 
    dataTbl[11, 1] <- glabel("Dur. years", container = dataTbl)
    dataTbl[11, 2] <- (datW$rMax1Dur <- gedit(" ", width = 50, container = dataTbl))
    
    dataTbl[12, 1:2] <- gseparator(horizontal = TRUE, container = dataTbl, expand = TRUE)
    
    ## rMax2 data 
    dataTbl[13, 1 ] <- glabel("Hist. Max 2", container = dataTbl)
    dataTbl[13, 2 ] <- (datW$rMax2 <- gedit(" ", width = 50, container = dataTbl))
    
    ## rMax2 duration 
    dataTbl[14, 1] <- glabel("Dur. years", container = dataTbl)
    dataTbl[14, 2] <- (datW$rMax2Dur <- gedit(" ", width = 50, container = dataTbl))
    
    ## valid all fields
    dataTbl[15, 2] <- (gd <- ggroup(container = dataTbl))
    (datW$OK <- gbutton("OK", container =gd, handler = bindAllData))
    
    datW$sb <- gstatusbar("choose a csv file and edit Project/Main if needed",
                          container = renGUI)
    
    ## updateData()
    addhandlerchanged(obj = datW$project,
                      handler = updateData, action = datW)
    
    addHandlerChanged(obj = datW$rMax1,
                      handler = validNumList, action = datW)
    
    addHandlerChanged(obj = datW$rMax2,
                      handler = validNumList, action = datW)
    
    ## Use the help line
    
    ## addHandlerMouseMotion(obj = datW[["New"]], handler = explain,
    ##                      action = "Create new project from csv file")
    
    ## 2011-12-13
    ## addHandlerMouseMotion(obj = dataTbl, handler = explain,
    ##                      action = "Edit the chosen project")
    
    addHandlerMouseMotion(obj = datW$main, handler = explain,
                          action = "Enter main title for data (only)")
    
    addHandlerMouseMotion(obj = datW$effDuration, handler = explain,
                          action = "Enter effective duration for the main sample (in years)")
    
    
    for (nm in c("rMax1Dur", "rMax2Dur")) {
        addHandlerMouseMotion(obj = datW[[nm]], handler = explain,
                              action = "Enter effective duration of the historical period (in years)")
    }
    
    for (nm in c("rMax1", "rMax2")) {
        addHandlerMouseMotion(obj = datW[[nm]], handler = explain,
                              action = "Enter maximal values over the period e.g. \"129.30;120.10;178.28\"")
    }
    
    ## for (objName in c("main", "effDuration", "rMax1", "rMax1Dur", "rMax2", "rMax2Dur")) {
    ##   addHandlerChanged(obj = datW[[objName]], handler = bindData, action = objName)
    ## }
    
    ## 2011-10-14 Add default values
    for ( name in names(re$RLprov$.defaultProj) ) {
        svalue(datW[[name]]) <- re$RLprov$.defaultProj[[name]]
        enabled(datW[[name]]) <- FALSE
    }
    
    enabled(datW[["OK"]]) <- FALSE
    
    addHandlerChanged(obj = datW[["main"]], handler = bindData, action = "main")
    addHandlerChanged(obj = datW[["effDuration"]], handler = bindData, action = "effDuration")
    addHandlerChanged(obj = datW[["rMax1"]], handler = bindData, action = "rMax1")
    addHandlerChanged(obj = datW[["rMax1Dur"]], handler = bindData, action = "rMax1Dur")
    addHandlerChanged(obj = datW[["rMax2"]], handler = bindData, action = "rMax2")
    addHandlerChanged(obj = datW[["rMax2Dur"]], handler = bindData, action = "rMax2Dur")
    
    
    if (DEBUG) cat("... done\n")
    
    ##===============================================================
    ## fit part
    ##=============================================================== 
    
    if (DEBUG) cat("Build 'fit' widgets...\n")
    
    supportDists <- c(exponential = "exp",
                      ## gpd = "gpd",
                      GPD = "GPD",
                      weibull = "weibull",
                      mixexp2 = "mixexp2")
    
    fitW <- list()
    
    fitTbl[1, 1:2] <- ( gfit0 <- ggroup(container =fitTbl) )
    
    gbutton(text = "New", container =gfit0, handler = newFit)
    gbutton(text = "Delete", container =gfit0, handler = delFit)
    
    fitTbl[2, 1] <- glabel("Fit", container =fitTbl)
    fitTbl[2, 2] <- (fitW$fit <- gcombobox(items = "           ",
                                           selected = 0, container = fitTbl))
    
    fitTbl[3, 1] <- glabel("Projet", container =fitTbl)
    fitTbl[3, 2] <- ( fitW$project <- glabel("         ", container = fitTbl) )
    
    fitTbl[4, 1] <- glabel("Eff. duration (yrs)", container =fitTbl)
    fitTbl[4, 2] <- ( fitW$effDuration <- glabel("", container =fitTbl,
                                                 coerce.with = as.numeric) )
    
    fitTbl[5, 1:3] <- gseparator(horizontal = TRUE, container = fitTbl, expand = TRUE)
    
    ##=================================================================================
    ## 'useMax1' and 'useMax2' are taken as gradios and not check boxes since
    ## enabled does not work for checkboses in the tcltk toolkit at the time.
    ##=================================================================================
    
    fitTbl[6, 1] <- glabel("Use Hist. 1", container =fitTbl) 
    fitTbl[6, 2] <- ( fitW$useMax1 <- gradio(items = c("yes", "no"), selected = 1,
                                             horizontal = TRUE, container =fitTbl) )
    
    fitTbl[6, 3] <- ( fitW$rMax1Lab <- glabel("                          ", container =fitTbl) )
    
    fitTbl[7, 1] <- glabel("Use Hist. 2", container =fitTbl) 
    fitTbl[7, 2] <- ( fitW$useMax2 <- gradio(items = c("yes", "no"), selected = 1,
                                             horizontal = TRUE, container =fitTbl) )
    
    fitTbl[7, 3] <- ( fitW$rMax2Lab <- glabel("                ", container =fitTbl) )
    
    enabled(fitW$useMax1) <- FALSE
    enabled(fitW$useMax2) <- FALSE
    enabled(fitW$rMax1Lab) <- FALSE
    enabled(fitW$rMax2Lab) <- FALSE
    
    ##---------------------------------------------------------------------------------
    fitTbl[8, 1:3] <- gseparator(horizontal = TRUE, container = fitTbl, expand = TRUE)
    
    ## Threshold part
    fitTbl[9, 1] <- glabel("Fit threshold", container =fitTbl)
    fitTbl[9, 2] <-  (fitW$threshold <- gspinbutton(from = 0, to = 1, by = 0.1, container =fitTbl) )
    
    ## Distribution part
    fitTbl[10, 1] <- glabel("Distr.", container =fitTbl)
    fitTbl[10, 2:3] <- ( fitW$distribution  <- gradio(items = names(supportDists),
                                                      selected = 1,
                                                      container = fitTbl,
                                                      horizontal = TRUE) )
    
    fitTbl[11, 1] <- glabel("Conf. lev.", container =fitTbl)
    fitTbl[11, 2:3] <- ( fitW$confLev <- gedit("95;70",
                                               width = 50,
                                               container =fitTbl) )
    
    fitTbl[12, 1] <- glabel("Return periods", container =fitTbl)
    fitTbl[12, 2:3] <- ( fitW$returnPeriods <- gedit("100;200;500;1000",
                                                     width = 50,
                                                     container =fitTbl) )
    
    ##fitTbl[13, 1] <- glabel("Main title", container =fitTbl)
    ##fitTbl[13, 2:3] <- ( fitW$main <- gedit("", container =fitTbl) )
    
    
    ##fitTbl[6, 1] <- glabel("conf. lev.", container =fitTbl)
    fitTbl[14, 4] <- (fitW$OK <- gbutton("FIT", container =fitTbl))
    enabled(fitW$OK) <- FALSE
    
    addHandlerMouseMotion(obj = fitTbl, handler = explain,
                          action = "edit the fit options")
    
    addHandlerMouseMotion(obj = fitW$threshold, handler = explain,
                          action = "threshold for POT")
    
    addHandlerMouseMotion(obj = fitW$distribution, handler = explain,
                          action = "distribution for POT exceedances")
    
    addHandlerMouseMotion(obj = fitW$confLev, handler = explain,
                          action = "percents for confidence intervals")
    
    addHandlerMouseMotion(obj = fitW$returnPeriods, handler = explain,
                          action = "return periods in years for results")
    
    ## addHandlerMouseMotion(obj = fitW$main, handler = explain,
    ##                      action = "main title for fits")
    
    addHandlerMouseMotion(obj = fitW$rMax1Lab, handler = explain,
                          action = "historical data 1 (not editable)")
    
    addHandlerMouseMotion(obj = fitW$rMax2Lab, handler = explain,
                          action = "historical data 2 (not editable)")
    ##
    addHandlerChanged(obj = fitW$useMax1, handler = bindFits, action = "useMax1")
    addHandlerChanged(obj = fitW$useMax2, handler = bindFits, action = "useMax2")
    addHandlerChanged(obj = fitW$threshold, handler = bindFits, action = "threshold")
    ## 2011-10-18
    addHandlerKeystroke(obj = fitW$threshold,
                        key = "e.ENTER",
                        handler = bindFits, action = "threshold")
    
    addHandlerChanged(obj = fitW$distribution, handler = bindFits, action = "distribution")
    
    ## for (objName in c("useMax1", "useMax2", "threshold", "distribution")) {
    ##   addHandlerChanged(obj = fitW[[objName]], handler = bindFits, action = objName)
    ## }
    
    ## DOES NOT WORK IN A LOOP WITH tcl-tk. Probably a question of evaluation/bindng in tcl
    addHandlerKeystroke(obj = fitW$confLev,
                        key = "e.ENTER",
                        handler = bindFits, action = "confLev")
    addHandlerKeystroke(obj = fitW$returnPeriods,
                        key = "e.ENTER",
                        handler = bindFits, action = "returnPeriods")
    
    
    ## 2011-10-14 Add default values
    for ( name in names(re$RLprov$.defaultFit) ) {
        ## cat("XXX setting name = ", name, "value = ", RLprov$.defaultFit[[name]], "\n") 
        svalue(fitW[[name]]) <- re$RLprov$.defaultFit[[name]]
        enabled(fitW[[name]]) <- FALSE
    }
    
    enabled(fitW[["OK"]]) <- FALSE
    
    
    ## ABANDONNED but efficient in RGktk2
    ## for (objName in c("returnPeriods", "main")) {
    ##   addHandlerKeystroke(obj = fitW[[objName]], key = "e.ENTER", handler = bindFits, action = objName)  
    ## }
    
    if (DEBUG) cat("... done\n")
    
    ##===============================================================
    ## results part
    ##=============================================================== 
    
    if (DEBUG) cat("Build 'results' widgets...\n")
    
    resW <- list()
    
    resTbl <- glayout(container = resGrp, label = "")
    ## resW[["NB"]]

    resTbl[1, 1] <- ( gr0 <- ggroup(container = resTbl, horizontal = TRUE) )
    
    ( resW$plot <- gbutton(text = "Plot",
                           active = FALSE,
                           handler = plotRes,
                           action = list(device = "x11"),
                           container = gr0) )
    
    ( resW$rlOptions <- gbutton(text = "Options",
                                handler = rlOptions,
                                container = gr0) )
    
    
    ( resW$delete <- gbutton(text = "Delete",
                             handler = delRes,
                             container = gr0) )
    
    ( resW$exportHTML <- gbutton(text = "Export HTML", container = gr0,
                                 handler = HTMLRes) )
    
    ## ( resW$exportODT <- gbutton(text = "Export ODT", container = gr0,
    ##                          handler = exportRes) )
    
    ( resW$dumpR <- gbutton(text = "Dump R", container = gr0,
                            handler = dumpRes) )
    
    resW$nB <- gnotebook(container = resGrp, expand = TRUE)
    
    ##
    if (FALSE) {
        
        addHandlerMouseMotion(obj = resW$exportODT, handler = explain,
                              action = "export results in odt (Open Office Text document) format. ")
        
        addHandlerMouseMotion(obj = resW$exportHTML, handler = explain,
                              action = "export results as HTML document with image link(s)")
        
        addHandlerMouseMotion(obj = resW$delete, handler = explain,
                              action = "delete active results")
        
        addHandlerMouseMotion(obj = resW$plot, handler = explain,
                              action = "plot return levels for active results")
        
        addHandlerMouseMotion(obj = resW$rlOptions, handler = explain,
                              action = "edit options for the current graph: title, ...")
        
        addHandlerMouseMotion(obj = resW$nB, handler = explain,
                              action = "active results")
    }
    
    Test <- data.frame(period = c(1.0), prob = c(0.0), quant = c(0.0),
                       L.95 = NA, U.95 = NA, L.70 = NA, U.70 = NA)
    
    ##===============================================================
    ## Add handler. This must be done only now due to scoping
    ## problems in widgets.
    ##=============================================================== 
    
    
    addhandlerchanged(obj = fitW$fit, handler = updateFit, action = fitW)
    
    if (DEBUG) cat("  'fitProj'...\n")
    addhandlerclicked(obj = fitW$OK, handler = fitProj, action = fitW)
    
    if (DEBUG) cat("  'updateRes'...\n")
    ## addhandlerchanged(obj = resW$fit, handler = updateRes)
    
    if (DEBUG) cat("... done\n")
    
    if (DEBUG) cat("Make window visible...\n")
    
    svalue(noteBook) <- 1
    
    size(renGUI) <- c(900, 600)
    visible(renGUI) <- TRUE
    
    if (DEBUG) {
        cat("main window size:", size(renGUI), "\n")
    }
    
    
    if (DEBUG) cat("... done\n")
    
    menuCopyCoef = list(copy = gaction(label = "copy coef",  handler = coef2Clipboard))
    menuCopyPred = list(copy = gaction(label = "copy pred", handler = pred2Clipboard))
    
}
