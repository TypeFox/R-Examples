library(shiny)
library(QCA)
# options(shiny.maxRequestSize=30*1024^2) # 30MB per file
# options(warn=-1)

setwd(Sys.getenv("userwd"))
options(width=74)
options(help_type = "html")
ev <- new.env()
## Future feature: local settings

# if (file.exists(noedit)) {
#     settings <- strsplit(readLines(noedit), split=" = ")
#     opts <- unlist(lapply(settings, "[[", 1))
#     vals <- unlist(lapply(settings, "[[", 2))
# }

listFiles <- function(dirpath, filetype="*") {
    
    result <- list(dirs = NULL, files = NULL, filepath = as.matrix(filepath), ok = oktoset)
    
    # get all files
    found <- list.files(dirpath)
    
    temp <- toupper(found)
    found <- found[match(sort(temp), temp)]
    
    if (length(found) > 0) {
        
        isdir <- sapply(found, function(x) file_test("-d", file.path(dirpath, x)))
        
        if (any(isdir)) {
            result$dirs <- found[isdir]
        }
        
        if (any(!isdir)) {
            if (filetype != "*") {
                extensions <- unlist(lapply(strsplit(found, split="\\."), "[", 2))
                found <- found[which(toupper(extensions) == toupper(filetype))]
                if (length(found) > 0) {
                    result$files <- found
                }
            }
            else {
                result$files <- found[!isdir]
            }
        }
        
        # this is necessary because Javascript interprets an array of length 1
        # as a singleton (not an array anymore)
        
        if (length(result$dirs) == 1) {
            result$dirs <- as.matrix(result$dirs)
        }
        
        if (length(result$files) == 1) {
            result$files <- as.matrix(result$files)
        }
        
    }
    
    # result$filename <- as.matrix(unlist(strsplit(basename(filepath), split="\\."))[1])
    if (!identical(filepath, "")) {
        if (file_test("-f", filepath)) {
            
            resfilename <- unlist(strsplit(basename(filepath), split="\\."))
            
            extension <<- resfilename[length(resfilename)]
            resfilename <- resfilename[-length(resfilename)]
            resfilename <- paste(gsub("[[:space:]]", "", gsub("[^[:alnum:] ]", "", resfilename)), collapse="")
            
            if (QCA::possibleNumeric(substr(resfilename, 1, 1))) {
                resfilename <- paste("x", resfilename, sep="")
            }
            
            filename <<- resfilename
            
        }
    }
    
    
    result$filename <- filename
    result$extension <- extension
    result$wd <- getwd()
    
    return(result)
    
}


tempdata <- NULL
mydata <- NULL
mytt <- NULL
myeqmcc <- NULL



current_path <- getwd()
oktoset <- TRUE
filepath <- ""
filename <- ""
extension <- ""
tcisdata <- TRUE

# vertical random noise for the thresholds setter in the calibrate dialog
vert <- NULL



# Define server logic for random distribution application
shinyServer(function(input, output, session) {
    
    observe({
        dirfilist <- input$dirfilist
        if (!is.null(dirfilist)) {
            if (dirfilist$refresh) {
                filename <<- ""
            }
        }
        session$sendCustomMessage(type = "dirfile", listFiles(current_path))
    })
    
    
    # import csv files
    observe({
        
        # this changes every time the Javascript object "read_table" changes
        read_table <- input$read_table
        
        filepath <<- ""
        
        oktoset <<- TRUE
            
        
        if (!is.null(input$dirfile_chosen)) {
            
            dfchosen <- input$dirfile_chosen
            oktoset <<- TRUE
            
            # run this each time the "dirfile_chosen" vector changes, see Javascript code
            
            if (dfchosen[1] == "file") {
                
                filepath <<- file.path(gsub("[/]$", "", current_path), dfchosen[2])
                
            }
            else {
                splitpath <- unlist(strsplit(current_path, split=.Platform$file.sep))
                
                if (dfchosen[2] == ".." | dfchosen[2] == "...") {
                    
                    if (length(splitpath) > 1) {
                        splitpath <- splitpath[-length(splitpath)]
                    }
                    
                    if (identical(splitpath, "")) {
                        splitpath <- "/"
                    }
                    
                    pathtobe <- paste(splitpath, collapse=.Platform$file.sep)
                    
                    # if (file_test("-x", pathtobe)) {
                    if (length(list.files(pathtobe)) > 0) {
                        current_path <<- pathtobe
                    }
                    else {
                        oktoset <<- FALSE
                    }
                }
                else {
                    current_dirs <- listFiles(current_path)$dirs
                    
                    if (dfchosen[2] %in% current_dirs) {
                        pathtobe <- file.path(current_path, dfchosen[2])
                        pathtobe <- paste(pathtobe, "/", sep="")
                        
                        # if (file_test("-x", pathtobe)) {
                        if (length(list.files(pathtobe)) > 0) {
                            current_path <<- pathtobe
                        }
                        else {
                            oktoset <<- FALSE
                        }
                    }
                    else if (dfchosen[2] != "__stdir__") {
                        
                        if (dfchosen[2] == "root") {
                            dfchosen[2] <- ""
                        }
                        
                        splitpath <- splitpath[seq(which(splitpath == dfchosen[2]))]
                        pathtobe <- ifelse(length(splitpath) == 1,
                                           ifelse(identical(splitpath, ""), "/", splitpath),
                                           paste(splitpath, collapse=.Platform$file.sep))
                        
                        # if (file_test("-x", pathtobe)) {
                        if (length(list.files(pathtobe)) > 0) {
                            current_path <<- pathtobe
                        }
                        else {
                            oktoset <<- FALSE
                        }
                    }
                }
            }
            
            if (oktoset) {
                # This is mainl for the Unix systems,
                # where the root is a backslash
                if (!grepl("/", current_path)) {
                    current_path <<- paste(current_path, "/", sep="")
                }
            }
            
            if (dfchosen[1] == "dir" & dfchosen[3] != "") {
                
                if (grepl("cannot change working directory", tryCatch(setwd(dfchosen[3]), error = function(e) e))) {
                    listtosend <- listFiles(current_path)
                    listtosend$filename <- "error!"
                    session$sendCustomMessage(type = "dirfile", listtosend)
                }
                else {
                    
                    if (dfchosen[3] %in% listFiles(current_path)$dirs) {
                        pathtobe <- file.path(current_path, dfchosen[3])
                        
                        # if (file_test("-x", pathtobe)) {
                        if (length(list.files(pathtobe)) > 0) {
                            current_path <<- pathtobe
                        }
                    }
                    else {
                        current_path <<- dfchosen[3]
                    }
                    
                    # This is mainly for the Windows system,
                    # where a directory path needs to be terminated with a backslash
                    if (!grepl("/", current_path)) {
                        current_path <<- paste(current_path, "/", sep="")
                    }
                    
                    setwd(current_path)
                    session$sendCustomMessage(type = "dirfile", listFiles(current_path))
                }
            }
            else {
                setwd(current_path)
                session$sendCustomMessage(type = "dirfile", listFiles(current_path))
            }
                            
        }
        else {
            current_path <<- getwd()
        }
        
        
        if (!identical(filepath, "")) {
            
            numevars <- ""
            
            header <- read_table$header
            colsep <- read_table$sep # comma separated by default
            row_names <- read_table$row_names
            decimal <- read_table$dec
            
            filename <<- unlist(strsplit(basename(filepath), split="\\."))
            filename <<- filename[-length(filename)]
            if (length(filename) > 1) {
                filename <<- paste(filename, collapse=".")
            }
            
            
            if (!identical(row_names, "")) { # this isn't a vector, just a name or a number, but just as a best practise
                if (QCA::possibleNumeric(row_names)) {
                    row_names <- as.numeric(row_names)
                }
                tc <- capture.output(tryCatch(read.table(filepath, header=header, ifelse(colsep == "tab", "\t", colsep),
                          row.names=row_names, as.is=TRUE, dec=decimal, nrows = 2), error = function(e) e, warning = function(w) w))
                   
            }
            else {
                tc <- capture.output(tryCatch(read.table(filepath, header=header, ifelse(colsep == "tab", "\t", colsep),
                          as.is=TRUE, dec=decimal, nrows = 2), error = function(e) e, warning = function(w) w))
            }
            
            if (any(grepl("subscript out of bounds", tc))) {
                mesaj <- paste("The data doesn't have ", row_names, " columns.", sep = "")
                session$sendCustomMessage(type = "tempdatainfo", list(ncols=1, nrows=1, colnames=mesaj, rownames="error!"))
                return(invisible())
            }
            else if (any(grepl("are not allowed", tc))) {
                mesaj <- paste("The row.names column has duplicated values.", sep = "")
                session$sendCustomMessage(type = "tempdatainfo", list(ncols=1, nrows=1, colnames=mesaj, rownames="error!"))
                return(invisible())
            }
            else if (any(grepl("data frame with 0 columns", tc))) {
                mesaj <- paste("The data has only 1 column.", sep = "")
                tc <- tryCatch(read.table(filepath, header=header, ifelse(colsep == "tab", "\t", colsep),
                           as.is=TRUE, dec=decimal, nrows = 2), error = function(e) e)
                session$sendCustomMessage(type = "tempdatainfo", list(ncols=2, nrows=2, colnames=c(colnames(tc), mesaj), rownames=""))
                return(invisible())
            }
            else if (any(grepl("attempt to select less than one element", tc))) {
                mesaj <- paste("The column \"", row_names, "\" was not found.", sep = "")
                session$sendCustomMessage(type = "tempdatainfo", list(ncols=1, nrows=1, colnames=mesaj, rownames="error!"))
                return(invisible())
            }
            
            
            tc <- tryCatch(read.table(filepath, header=header, ifelse(colsep == "tab", "\t", colsep),
                           as.is=TRUE, dec=decimal, nrows = 2), error = function(e) e, warning = function(w) w)
            
            
            tcisdata <<- TRUE
            
            
            if (is.null(dim(tc))) {
                if (is.list(tc)) {
                    if (identical(names(tc), c("message", "call"))) {
                        tcisdata <<- FALSE
                        session$sendCustomMessage(type = "tempdatainfo", list(ncols=1, nrows=1, colnames=tc$message, rownames="error!"))
                    }
                }
            }
            else {
                if (grepl("X.PDF", names(tc)[1])) {
                    tcisdata <<- FALSE
                    session$sendCustomMessage(type = "tempdatainfo", list(ncols=1, nrows=1, colnames="not a dataframe, this is a PDF file", rownames="error!"))
                }
            }
            
            if (tcisdata) {
                
                if (row_names != "") {
                    tc <- tryCatch(read.table(filepath, header=header, ifelse(colsep == "tab", "\t", colsep),
                              row.names=row_names, as.is=TRUE, dec=decimal), error = function(e) e, warning = function(w) w)
                }
                else {
                    tc <- tryCatch(read.table(filepath, header=header, ifelse(colsep == "tab", "\t", colsep),
                              as.is=TRUE, dec=decimal), error = function(e) e, warning = function(w) w)
                }
                
                
                if (identical(names(tc), c("message", "call"))) {
                    # not a mistake, it might happen that a warning is issued even at this stage
                    session$sendCustomMessage(type = "tempdatainfo", list(ncols=1, nrows=1, colnames=tc$message, rownames="error!"))
                }
                else {
                
                    tempdata <<- tc
                    assign(filename, tc, envir = ev)
                    
                    session$sendCustomMessage(type = "tempdatainfo", list(ncols=ncol(tempdata),
                                                                      nrows=nrow(tempdata),
                                                                      colnames=colnames(tempdata),
                                                                      rownames=rownames(tempdata)))
                    
                }
            }
        }
        
    })
    
    
    observe({
        
        import <- input$import
        
        if (!is.null(import) & tcisdata) {
            
            mydata <<- tempdata
            
            set.seed(12345)
            vert <<- sample(185:200, nrow(mydata), replace = TRUE)
            
            numerics <- as.vector(unlist(lapply(mydata, QCA::possibleNumeric)))
            
            calibrated <- as.vector(unlist(lapply(mydata, function(x) {
                all(na.omit(x) >= 0 & na.omit(x) <= 1)
            })))
            
            rowend <- min(17, nrow(mydata))
            colend <- min(8, ncol(mydata))
            
            
            tosend <- as.list(mydata[seq(rowend), seq(colend)])
            names(tosend) <- NULL # to make it look like a simple Array in Javascript
            
            session$sendCustomMessage(type = "datainfo",
                                      list(list(ncols=ncol(mydata),              # the datainfo
                                                nrows=nrow(mydata),
                                                colnames=colnames(mydata),
                                                rownames=rownames(mydata),
                                                numerics=numerics,
                                                calibrated=calibrated),
                                            tosend,                               # theData
                                            paste(1, 1, rowend, colend, ncol(mydata), sep="_"))) # and the dataCoords
            
        }
    })
    
    
    
    observe({
        scrollvh <- input$scrollvh
        
        if (!is.null(scrollvh)) {
            
            scrollvh <- scrollvh + 1
            
            rowstart <- scrollvh[1]
            colstart <- scrollvh[2]
            rowend <- min(rowstart + scrollvh[3] - 1, nrow(mydata))
            colend <- min(colstart + scrollvh[4] - 1, ncol(mydata))
            
            tosend <- as.list(mydata[seq(rowstart, rowend), seq(colstart, colend)])
            names(tosend) <- NULL
            
            # session$sendCustomMessage(type = "theData", tosend)
            session$sendCustomMessage(type = "theData", list(tosend, paste(rowstart, colstart, rowend, colend, ncol(mydata), sep="_")))
        }   
    })
    
    
    
    
    # eqmcc
    observe({
        
        foo <- input$eqmcc2R
        
        if (!is.null(mydata)) {
            
            outc <- c("")
            if (length(foo$outcome) > 0) {
                outc <- unlist(foo$outcome)
            }
            
            cnds <- c("")
            if (length(foo$conditions) > 0) {
                cnds <- unlist(foo$conditions)
            }
            
            expl <- c("1")
            if (length(foo$explain) > 0) {
                expl <- unlist(foo$explain)
            }
            
            incl <- c("")
            if (length(foo$include) > 0) {
                incl <- unlist(foo$include)
            }
            
            direxp <- ""
            if (length(foo$dir_exp) > 0) {
                direxp <- unlist(foo$dir_exp)
                if (all(direxp == "-")) {
                    direxp <- ""
                }
            }
            
            use_letters <- foo$use_letters
            
            myeqmcc <<- NULL
            incl.cut <- c(as.numeric(foo$ic1), as.numeric(foo$ic0))
            
            textoutput <- capture.output(tryCatch(
                (myeqmcc <<- eqmcc(mydata, outcome = outc,
                      neg.out = foo$neg_out,
                      conditions = cnds,
                      relation = foo$relation,
                      n.cut = as.numeric(foo$n_cut),
                      incl.cut = incl.cut[!is.na(incl.cut)],
                      explain = expl,
                      include = incl,
                      all.sol = foo$all_sol,
                      dir.exp = direxp,
                      details = foo$details,
                      show.cases = foo$show_cases,
                      use.tilde = foo$use_tilde,
                      use.letters = use_letters,
                      via.web=TRUE)) , error = function(e) e)
                      # PRI=foo$PRI)) , error = function(e) e)
            )
            
            if (any(error <- grepl("Error", textoutput))) {
                errmessage <- paste0("Error:", unlist(strsplit(textoutput[which(error)], split=":"))[2])
                errmessage <- substr(errmessage, 1, nchar(errmessage) - 1)
                textoutput <- c("error", errmessage, "")
            }
            
            textoutput <- gsub("undefined columns selected>", "Column names in the command don't match those in the interface.", textoutput)
            
            sendnormal <- FALSE
            
            if (length(cnds) <= 7) { # to draw the Venn diagram, up to 5 variables for now
                
                if (!identical(outc, "")) { # this happens when the interface is disconnected
                    if (length(QCA::splitstr(outc)) == 1 & !is.null(myeqmcc)) {
                        
                        myeqmcc$tt$initial.data <- NULL
                        myeqmcc$tt$recoded.data <- NULL
                        myeqmcc$tt$indexes <- myeqmcc$tt$indexes - 1 # to take the indexes in Javascript notation
                        if (identical(cnds, "")) {
                            cnds <- toupper(setdiff(names(mydata), outc))
                        }
                        
                        if (use_letters) {
                            cnds <- LETTERS[seq(length(cnds))]
                        }
                        
                        myeqmcc$tt$options$conditions <- toupper(cnds)
                        myeqmcc$tt$id <- apply(myeqmcc$tt$tt[, toupper(cnds)], 1, function(x) {
                            ifelse(any(x == 1), paste(which(x == 1), collapse=""), "0")
                        })
                        
                        session$sendCustomMessage(type = "eqmcc", list(textoutput, list(as.list(myeqmcc$tt      ))))
                    }
                    else {
                        sendnormal <- TRUE
                    }
                }
                else {
                    sendnormal <- TRUE
                }
            }
            
            if (sendnormal) {
                session$sendCustomMessage(type = "eqmcc", list(textoutput, NULL))
            }
            
        }
        
    })
    
    
    # truth table
    observe({
        
        foo <- input$tt2R
        
        if (!is.null(mydata)) {
            
            outc <- ""
            if (length(foo$outcome) > 0) {
                outc <- unlist(foo$outcome)
            }
            
            cnds <- ""
            if (length(foo$conditions) > 0) {
                cnds <- QCA::splitstr(unlist(foo$conditions))
            }
            
            sortbys <- unlist(foo$sort_by)
            selected <- unlist(foo$sort_sel)
            sortbys <- sortbys[selected[names(sortbys)]]
            
            if (length(sortbys) == 0) {
                sortbys <- ""
            }
            
            use_letters <- foo$use_letters
            
            mytt <<- NULL
            
            incl.cut <- c(as.numeric(foo$ic1), as.numeric(foo$ic0))
            
            textoutput <- capture.output(tryCatch(
                (mytt <<- truthTable(mydata, outcome = outc,
                      neg.out = foo$neg_out,
                      conditions = cnds,
                      n.cut = as.numeric(foo$n_cut),
                      incl.cut = incl.cut[!is.na(incl.cut)],
                      complete = foo$complete,
                      show.cases = foo$show_cases,
                      sort.by = sortbys,
                      use.letters = foo$use_letters)), error = function(e) e)
                      # PRI=foo$PRI)), error = function(e) e)
            )
            
            if (any(error <- grepl("Error", textoutput))) {
                errmessage <- paste0("Error:", unlist(strsplit(textoutput[which(error)], split=":"))[2])
                errmessage <- substr(errmessage, 1, nchar(errmessage) - 1)
                textoutput <- c("error", errmessage, "")
            }
            
            textoutput <- gsub("undefined columns selected>", "Column names in the command don't match those in the interface.", textoutput)
            
            if (length(cnds) <= 7) { # to draw the Venn diagram, up to 5 variables for now
                if (!is.null(mytt)) {
                    mytt$initial.data <- NULL
                    mytt$recoded.data <- NULL
                    mytt$indexes <- mytt$indexes - 1 # to take the indexes in Javascript notation
                    if (identical(cnds, "")) {
                        cnds <- toupper(setdiff(names(mydata), outc))
                    }
                    
                    if (use_letters) {
                        cnds <- LETTERS[seq(length(cnds))]
                    }
                    
                    mytt$options$conditions <- toupper(cnds)
                    mytt$id <- apply(mytt$tt[, toupper(cnds)], 1, function(x) {
                        ifelse(any(x == 1), paste(which(x == 1), collapse=""), "0")
                    })
                    
                }
                
                session$sendCustomMessage(type = "tt", list(textoutput, mytt))
            }
            else {
                session$sendCustomMessage(type = "tt", list(textoutput, NULL))
            }
            
            if (!is.null(mytt)) {
                assign("tt", mytt, envir = ev)
            }
        }
    })
    
    
    
    observe({
        thinfo <- input$thinfo
        
        if (!is.null(thinfo)) {
            if (thinfo[2] != "") {
                
                # make sure the variable is numeric
                if (QCA::possibleNumeric((mydata[, thinfo[2]]))) {
                    response <- findTh(mydata[, thinfo[2]], n = as.numeric(thinfo[1]))
                }
                else {
                    response <- "notnumeric"
                }
                
                session$sendCustomMessage(type = "thvalsfromR", 
                    as.list(response)
                )
            }
        }
    })
    
    
    # export to file
    observe({
        exportobj <- input$exportobj
        
        if (!is.null(exportobj)) {
            
            filesep <- exportobj$sep
            if (filesep == "tab") {
                filesep <- "\t"
            }
            
            separator <- exportobj$sep
            filetowrite <- file.path(current_path, exportobj$filename)
            
            if (exportobj$newfile) {
                if (exportobj$filename != "") {
                    filetowrite <- file.path(current_path, exportobj$filename)
                }
            }
            
            export(mydata, filetowrite, sep=filesep, col.names=exportobj$header, caseid=exportobj$caseid)
        }
    })
    
    
    
    # calibrate
    observe({
        foo <- input$calibrate
        
        if (!is.null(foo)) {
            
            checks <- rep(TRUE, 9)
            
            checks[1] <- !is.null(mydata)
            # checks[2] <- length(foo$thresholds) > 0
            
            foo$thresholds <- unlist(foo$thresholds)
            
            
            nms <- unlist(foo$thnames)[foo$thresholds != ""]
            foo$thresholds <- foo$thresholds[foo$thresholds != ""]
            
            
            thrs <- suppressWarnings(as.numeric(foo$thresholds))
            if (any(!is.na(thrs))) {
                foo$thresholds <- as.numeric(thrs[!is.na(thrs)])
            }
            
            # if (checks[2]) {
            #     checks[3] <- all(!is.na(suppressWarnings(as.numeric(foo$thresholds))))
            # }
            
            # if (checks[3]) {
            #     foo$thresholds <- as.numeric(foo$thresholds)
            # }
            
            if (all(foo$thresholds == "")) {
                foo$thresholds <- NA
            }
            else {
                if (!is.null(nms)) {
                    if (length(nms) == length(foo$thresholds)) {
                        names(foo$thresholds) <- nms
                    }
                    else {
                        foo$thresholds <- NA
                    }
                }
            }
            
            foo$x <- unlist(foo$x)
            
            checks[4] <- !is.null(foo$x)
            if (checks[4]) {
                checks[5] <- foo$x != ""
            }
            
            
            if (checks[1] & checks[5]) {
                
                if (foo$x %in% names(mydata)) {
                    checks[6] <- is.numeric(mydata[, foo$x])
                }
            }
            
            checks[7] <- QCA::possibleNumeric(foo$idm)
            if (checks[7]) {
                foo$idm <- as.numeric(foo$idm)
            }
            
            checks[8] <- QCA::possibleNumeric(foo$below)
            if (checks[8]) {
                foo$below <- as.numeric(foo$below)
            }
            
            checks[9] <- QCA::possibleNumeric(foo$above)
            if (checks[9]) {
                foo$above <- as.numeric(foo$above)
            }
            
            scrollvh <- unlist(foo$scrollvh)
            scrollvh <- scrollvh + 1
            
            if (all(checks)) {
                
                textoutput <- capture.output(tryCatch(
                    calibrate(
                        mydata[, foo$x],
                        type = foo$type,
                        thresholds = foo$thresholds,
                        include = foo$include,
                        logistic = foo$logistic,
                        idm = foo$idm,
                        ecdf = foo$ecdf,
                        below = foo$below,
                        above = foo$above), error = function(e) e)
                )
                
                response <- vector(mode="list", length = 3)
                
                if (any(error <- grepl("Error", textoutput))) {
                    errmessage <- paste0("Error:", unlist(strsplit(textoutput[which(error)], split=":"))[2])
                    errmessage <- substr(errmessage, 1, nchar(errmessage) - 1)
                    textoutput <- c("error", errmessage, "")
                }
                else {
                    textoutput <- "no problem"
                    
                    mydata[, ifelse(foo$newvar != "", foo$newvar, foo$x)] <<- calibrate(
                            mydata[, foo$x],
                            type = foo$type,
                            thresholds = foo$thresholds,
                            include = foo$include,
                            logistic = foo$logistic,
                            idm = foo$idm,
                            ecdf = foo$ecdf,
                            below = foo$below,
                            above = foo$above)
                    
                    rowstart <- scrollvh[1]
                    colstart <- scrollvh[2]
                    rowend <- min(rowstart + scrollvh[3] - 1, nrow(mydata))
                    colend <- min(colstart + scrollvh[4] - 1, ncol(mydata))
                    
                    tosend <- as.list(mydata[seq(rowstart, rowend), seq(colstart, colend)])
                    names(tosend) <- NULL
                    
                    numerics <- as.vector(unlist(lapply(mydata, QCA::possibleNumeric)))
                    response[[2]] <- list(ncols=ncol(mydata),
                                          nrows=nrow(mydata),
                                          colnames=colnames(mydata),
                                          rownames=rownames(mydata),
                                          numerics=numerics)
                    response[[3]] <- list(tosend, paste(rowstart, colstart, rowend, colend, ncol(mydata), sep="_"))
                    
                    
                }
                
                response[[1]] <- as.list(textoutput)
                
                session$sendCustomMessage(type = "calibrate", response)
                
            }
        }
        
    })
    
    
    
    
    
    # recode
    observe({
        foo <- input$recode
        
        if (!is.null(foo)) {
            
            checks <- rep(TRUE, 3)
            
            checks[1] <- !is.null(mydata)
            
            foo$x <- unlist(foo$x)
            
            checks[2] <- !is.null(foo$x)
            
            if (checks[2]) {
                checks[3] <- foo$x != ""
            }
            
            
            foo$oldv <- unlist(foo$oldv)
            foo$newv <- unlist(foo$newv)
            uniques <- unique(foo$newv)
            
            rules <- ""
            for (i in seq(length(uniques))) {
                part <- paste(paste(foo$oldv[foo$newv == uniques[i]], collapse = ","), uniques[i], sep="=")
                rules <- paste(rules, part, ifelse(i == length(uniques), "", "; "), sep="")
            }
            
            
            scrollvh <- unlist(foo$scrollvh)
            scrollvh <- scrollvh + 1
            
            if (all(checks)) {
                
                textoutput <- capture.output(tryCatch(
                    recode(mydata[, foo$x],
                           rules = rules), error = function(e) e)
                )
                
                response <- vector(mode="list", length = 3)
                
                if (any(error <- grepl("Error", textoutput))) {
                    errmessage <- paste0("Error:", unlist(strsplit(textoutput[which(error)], split=":"))[2])
                    errmessage <- substr(errmessage, 1, nchar(errmessage) - 1)
                    textoutput <- c("error", errmessage, "")
                }
                else {
                    textoutput <- "no problem"
                    
                    mydata[, ifelse(foo$newvar != "", foo$newvar, foo$x)] <<- recode(
                            mydata[, foo$x],
                            rules = rules)
                    
                    rowstart <- scrollvh[1]
                    colstart <- scrollvh[2]
                    rowend <- min(rowstart + scrollvh[3] - 1, nrow(mydata))
                    colend <- min(colstart + scrollvh[4] - 1, ncol(mydata))
                    
                    tosend <- as.list(mydata[seq(rowstart, rowend), seq(colstart, colend)])
                    names(tosend) <- NULL
                    
                    numerics <- as.vector(unlist(lapply(mydata, QCA::possibleNumeric)))
                    response[[2]] <- list(ncols=ncol(mydata),
                                          nrows=nrow(mydata),
                                          colnames=colnames(mydata),
                                          rownames=rownames(mydata),
                                          numerics=numerics)
                    response[[3]] <- list(tosend, paste(rowstart, colstart, rowend, colend, ncol(mydata), sep="_"))
                    
                    
                }
                
                response[[1]] <- as.list(textoutput)
                
                session$sendCustomMessage(type = "recode", response)
                
            }
        }
        
    })
    
    
    
    # dataModif
    observe({
        val <- input$dataModif
        if (!is.null(val)) {
            if (val[1] == "r") {
                rownames(mydata)[as.numeric(val[2])] <<- val[3]
            }
            else if (val[1] == "c") {
                colnames(mydata)[as.numeric(val[2])] <<- val[3]
            }
            else {
                valn <- suppressWarnings(as.numeric(val[3]))
                mydata[as.numeric(val[1]), as.numeric(val[2])] <<- ifelse(!is.na(valn), valn, val[3])
            }
        }
    })
    
    
    
    # dataPoints
    observe({
        foo <- input$thsetter2R
        if (!is.null(foo)) {
            
            cond <- unlist(foo$cond)
            
            if (cond %in% names(mydata)) {
                horz <- sort(mydata[, cond])
                response <- cbind(horz[!is.na(horz)], vert[!is.na(horz)])
                session$sendCustomMessage(type = "dataPoints", response)
            }
        }
    })
    
    
    
    # XYplotPoints
    observe({
        foo <- input$xyplot
        
        if (!is.null(foo)) {
            
            if (all(c(foo$x, foo$y) %in% names(mydata))) {
                
                X <- mydata[, foo$x]
                Y <- mydata[, foo$y]
                
                rpofsuf <- list(pof(    X,     Y, rel = "suf"),
                                pof(1 - X,     Y, rel = "suf"),
                                pof(    X, 1 - Y, rel = "suf"),
                                pof(1 - X, 1 - Y, rel = "suf"))
                
                rpofsuf <- lapply(rpofsuf, function(x) {
                    frmted <- formatC(c(x$incl.cov$incl, x$incl.cov$cov.r, x$incl.cov$PRI), format="f", digits = 3)
                    return(frmted)
                })
                
                rpofnec <- list(pof(    X,     Y, ron = TRUE),
                                pof(1 - X,     Y, ron = TRUE),
                                pof(    X, 1 - Y, ron = TRUE),
                                pof(1 - X, 1 - Y, ron = TRUE))
                
                rpofnec <- lapply(rpofnec, function(x) {
                    frmted <- formatC(c(x$incl.cov$incl, x$incl.cov$cov.r, x$incl.cov$RoN), format="f", digits = 3)
                    return(frmted)
                })
                
                
                response = list(rownames(mydata),
                                mydata[, foo$x],
                                mydata[, foo$y],
                                rpofsuf,
                                rpofnec)
                session$sendCustomMessage(type = "xyplot", response)
                
            }
        }
    })
    
    
    observe({
        
        foo <- input$Rcommand
        
        if (!is.null(foo)) {
        
            tc <- tryCatch(eval(parse(text = foo[2]), envir = ev), error = function(e) e, warning = function(w) w)
            
            if (inherits(tc, "error")) {
                tc <- unlist(strsplit(as.character(tc), split = ": "))
                tosend <- list("error", gsub("\n", "", tc[2]), "")
            }
            else if (inherits(tc, "warning")) {
                tc <- unlist(strsplit(as.character(tc), split = ": "))
                co <- capture.output(tryCatch(eval(parse(text = foo[2]), envir = ev)))
                
                tosend <- c("warning", strsplit(co[1], split = ","), gsub("\n", "", tc[2]))
            }
            else {
                co <- capture.output(tryCatch(eval(parse(text = foo[2]), envir = ev)))
                tosend <- strsplit(co, split = ",")
            }
            
            session$sendCustomMessage(type = "Rcommand", tosend)
        }
        
        
        
        
        
        
    })
    
    
    
    
    observe({
        
        foo <- input$changes
        
        if (!is.null(foo)) {
        
            session$sendCustomMessage(type = "getChanges", readLines(system.file("ChangeLog", package="QCAGUI")))
        
        }
        
    })
    
    
    
    
    observe({
        
        foo <- input$help
        
        if (!is.null(foo)) {
        
            # startDynamicHelp(start = TRUE)
            browseURL(file.path(path.package("QCAGUI"), "staticdocs", "index.html"))
            # help("QCAGUI")
            
        }
        
    })
})

