##==============================================================================
## DO NOT 'ROXYGEN' THIS
## Most functions are private (non-exported) functions
##============================================================================== 


##=============================================================================
##' For a vector of datetime \code{x} and a set of time
##' non-overlapping periods, find if the elements of \code{x} fall in
##' one period.
##'
##' @title Indicator for the union of periods
##' 
##' @param x \code{POSIXct}.
##' 
##' @param periods Data frame of periods with two \code{POSIXct} columns
##' \code{start} and \end{} and one row by period.
##' 
##' @param check Passed to \code{timeints2bounds} for control.
##'
##' @return Logical vector with the same length as \code{x}.
##'
##' @author yves
##' 
##' @note The periods must be clean: non-overlapping and in ascending
##' order. The intervals are assumed to be left closed [ (  as in the
##' 'findInterval' function.
##' 
inInts <- function(x,
                   periods,
                   check = TRUE) {
  
    chgs <- timeints2bounds(periods[ , c("start", "end")], check = check)  
    period <- findInterval(x = as.numeric(x),
                           vec = as.numeric(chgs),
                           rightmost.closed = TRUE,
                           all.inside = FALSE)
    
    inInt <-  as.logical(period %% 2)
    
}

##==============================================================================
##' Build time intervals thar are non-overlapping and in ascending
##' order from a set of periods which may not fulfill these two
##' conditions.
##' 
##' @title Build non-overlapping time intervals in ascending order
##'
##' @param data A data frame with \code{POSIXct} columns 'start and 'end',
##' and a character colum 'comment'.
##'
##' @param start Class \code{POSIXct} or coerced to this class with
##' \code{tz = TZ}.
##'
##' @param end Class \code{POSIXct} or coerced to this class with
##' \code{tz = TZ}.
##'
##' @param TZ The time zone.
##'
##' @param trace Level of verbosity.
##'
##' @return A data frame 
##'
##' @author yves
##'
##' @note Changed on 2015-04-18: remind that the concatenation methods 'c' drops
##' the timezone attribute of POSIXct objects!
##' 
cleanInt <- function(data,
                     start = NULL,
                     end = NULL,
                     TZ = "GMT",
                     trace = 0) {
    
    ## limited control
    if (!("start" %in% colnames(data)) || !("end" %in% colnames(data))) {
        stop("'data' must contain columns 'start' and 'end'")
    }
    
    if (!(is(data$start, "POSIXct")) || !(is(data$end, "POSIXct"))) {
        stop("columns 'start' and 'end' must be of class \"POSIXct\"")
    }
    
    ind <- (data$start > data$end)
    if ( any(ind) ) {
        warning("'data' contains rows with 'start' > 'end'")
        data <- data[!ind, ]
    } 
    
    inds <- order(data$start)
    
    if (nrow(data) > 1 ) {
        
        if (trace) cat("simplifying periods...\n")
        
        startCol <- data[inds, ]$start
        endCol <- data[inds, ]$end
        commentCol <- data[inds, ]$comment
        
        start.old <- startCol[1]
        end.old <- endCol[1]
        comment.old <- commentCol[1]
        
        n <- 0
        
        for (i in 2L:length(startCol)){
            
            if (startCol[i] > end.old) {
                if (n == 0) {
                    ## start.clean <- start.old
                    ## end.clean <- end.old
                    ## comment.clean <- comment.old
                    daf <- data.frame(start = start.old,
                                      end = end.old,
                                      comment = comment.old,
                                      stringsAsFactors = FALSE)
                } else {
                    ## start.clean <- c(start.clean, start.old)
                    ## end.clean <- c(end.clean, end.old)
                    ## comment.clean <- c(comment.clean, comment.old)
                    daf <- rbind(daf,
                                 data.frame(start = start.old,
                                            end = end.old,
                                            comment = comment.old,
                                            stringsAsFactors = FALSE))
                }
                n <- n + 1
                start.old <- startCol[i]
                end.old <- endCol[i]
                comment.old <- commentCol[i]
            } else {
                end.old <- endCol[i]
            }
            
        }
        
        ## build a dataframe. It can not have zero rows at this stage,
        ## but care is needed if 'start' or 'end' formal is given
        
        ## start.clean <- c(start.clean, start.old)
        ## end.clean <- c(end.clean, end.old)
        ## comment.clean <- c(comment.clean, comment.old)
        daf <- rbind(daf, data.frame(start = start.old,
                                     end = end.old,
                                     comment = comment.old,
                                     stringsAsFactors = FALSE))
        
        ## daf <- data.frame(start = as.POSIXct(start.clean),
        ##                   end = as.POSIXct(end.clean),
        ##                  comment = I(comment.clean))
        
        if (as.numeric(trace) >= 1) {
            print(daf)
        }
        
    } else {
        daf <- data
    }
    
    ## Cut at begining or end if necessary
    if (!is.null(start)) {
        start <- as.POSIXct(start, tz = TZ)
        
        if (start > daf$start[1]) {
            ## period is a period number
            chgs <- timeints2bounds(daf[ , c("start", "end")])
            period <- findInterval(x = as.numeric(start),
                                   vec = as.numeric(chgs),
                                   rightmost.closed = TRUE,
                                   all.inside = FALSE)
            
            ## start is in a gap... or not
            inGap <- period %% 2
            numGap <- period %/% 2 + 1
            if (trace) {
                cat("restrict with 'start' = ", format(start), "inGap = ",
                    inGap, "numGap = ", numGap, "\n")
            }
            daf <- daf[ 1:nrow(daf) >= numGap, ]
            if ( (nrow(daf) > 0) && inGap) daf$start[1] <- start
        }
    }
    
    if (as.numeric(trace) >= 1) {
        print(daf)
    }
    
    if (nrow(daf) && (!is.null(end))) {
        end <- as.POSIXct(end, tz = TZ)
        if (end < daf$end[nrow(daf)]) {
            ## period is a period number
            chgs <- timeints2bounds(daf[ , c("start", "end")])
            period <- findInterval(x = as.numeric(end),
                                   vec = as.numeric(chgs),
                                   rightmost.closed = TRUE,
                                   all.inside = FALSE)
            ## end is in a gap... or not
            inGap <- period %% 2
            numGap <- period %/% 2 + 1 + inGap
            if (trace) {
                cat("restrict with 'end' = ", format(end), "inGap = ",
                    inGap, "numGap = ", numGap, "\n")
            }
            daf <- daf[1L:nrow(daf) < numGap, ]
            if ((nrow(daf) > 0) && inGap) daf$end[nrow(daf)] <- end
        }
    }
    
    daf
    
}

##==============================================================================
##' Parse a node with \code{event} children in a XML file.
##' 
##' @title Parse a node with \code{event} children
##'
##' @param node An XML node with \code{event} children.
##'
##' @param varName Name of the variable.
##'
##' @param TZ Time zone.
##'
##' @param trace Level of verbosity.
##'
##' @return \code{NULL} if the provided node has no \code{event}
##' child, and else a data frame with one \code{POSIXct} column
##' \code{date} a \code{Var} column for the variable, and a character
##' column \code{comment}.
##'
##' @author yves
##'
##' @note \code{date} is a \code{POSIXct} object, The format given in
##' the XML file is assumed to be understandable by R. Other formats
##' (such as locale specific) should if needed be used within csv
##' files and not within the XML file directly.
##' 
parse.eventsNode <- function(node,
                             varName = "Var",
                             TZ = "GMT",
                             trace = 1) {
    
    event.nodes   <- XML::xmlElementsByTagName(node, "event")
    if (length(event.nodes)>0) {  
        date <- NULL
        Var <- NULL
        comment  <- NULL
        
        ## 'date' is kept as character here
        for (i in 1L:length(event.nodes)) {
            node.i <- event.nodes[[i]]
            date <- c(date, XML::xmlGetAttr(node.i, "date"))
            Var <- c(Var, as.numeric(XML::xmlValue(node.i)))
            comment.i <- XML::xmlGetAttr(node.i, "comment")
            if (length(comment.i) == 0) comment.i <- ""
            comment <- c(comment, comment.i)
        }
        
        ## 'date' is kept as character here. Are there true
        ## datetime in it?
        date2 <- as.POSIXct(rep(NA, length(date)), tz = TZ)
        ind <- grep(paste("[0-9]{1,4}-[0-1][0-9]-[0-3][0-9]($|( ", 
                          "([0-1][0-9]|[2][0-3]):([0-5][0-9])))", sep = ""),
                    date, perl = TRUE)
       
        if (length(ind)) {
            date[!ind] <- NA
            date2 <- as.POSIXct(date, tz = TZ)
        } else {
            date2 <- as.POSIXct(rep(NA, length(date)), tz = TZ)
        }
        
        Data <- data.frame(date = date2,
                           Var = as.numeric(Var),
                           comment = comment,
                           stringsAsFactors = FALSE)
        colnames(Data)[2] <- varName
        Data
    } else {
        NULL
    }
    
}

##==============================================================================
##' Parse a node with \code{period} children in a XML file.
##' 
##' @title Parse a node with \code{period} children
##'
##' @param node A XML node with \code{period} children.
##'
##' @param TZ Time zone.
##'
##' @param trace Level of verbosity.
##'
##' @return \code{NULL} if the provided node has no \code{period}
##' child, and else a data frame with two \code{POSIXct} columns
##' \code{start} and \code{end}, and a character column
##' \code{comment}.
##'
##' @author yves
##'
##' @note \code{date} is a \code{POSIXct} object, The format given in
##' the XML file is assumed to be understandable by R. Other formats
##' (such as locale specific) should if needed be used within csv
##' files and not within the XML file directly.
##' 
parse.periodsNode <- function(node,
                              TZ = "GMT",
                              trace = 1) {
    
    stopifnot(requireNamespace("XML", quietly = TRUE))
    
    period.nodes   <- XML::xmlElementsByTagName(node, "period")
    ## year.nodes   <- XML::xmlElementsByTagName(node, "year")
    if (length(period.nodes)>0) {
        start <- NULL
        end <- NULL
        comment  <- NULL
        
        for (i in 1L:length(period.nodes)) {
            node.i <- period.nodes[[i]]
            start <- c(start, XML::xmlGetAttr(node.i, "start"))
            end <- c(end, XML::xmlGetAttr(node.i, "end"))
            comment.i <- XML::xmlGetAttr(node.i, "comment")
            if (length(comment.i) == 0) comment.i <- ""
            comment <- c(comment, comment.i)
        }
       
        Prov <- data.frame(start = as.POSIXct(start, tz = TZ),
                           end = as.POSIXct(end, tz = TZ),
                           comment = as.character(comment),
                           stringsAsFactors = FALSE)
        
        Data <- cleanInt(data = Prov, trace = trace)
        
    } else{
        NULL
    } 
    
}

##=============================================================================
##' Read an 'events' file with description given in \code{node}.
##'
##' The varname must be given since it is an attribute of an ancestor
##' of the given node.
##' 
##' @title Read an 'events' file with description given in \code{node}
##'
##' @param node An XML node with attributes telling the number of columns,
##' which is the variable, and so on.
##'
##' @param dir Path of the directory to read from.
##'
##' @param varName The variable name.
##'
##' @param TZ Time zone.
##'
##' @param trace 
##'
##' @return A data frame with the events read and suitable columns.
##'
##' @author yves
##' 
read.eventsFile <- function(node,
                            dir,
                            varName = "Var",
                            TZ = "GMT",
                            trace = 1) {

    stopifnot(requireNamespace("XML", quietly = TRUE)) 
    path <- file.path(dir, XML::xmlGetAttr(node, "path"))
    
    if (trace) cat("Reading the data file\n", path,"...\n")
    nbCols <- as.integer(XML::xmlGetAttr(node, "nbCols"))
    dtCol <- as.integer(XML::xmlGetAttr(node, "dtCol"))
    varCol <- as.integer(XML::xmlGetAttr(node, "varCol"))
    commentCol <- as.integer(XML::xmlGetAttr(node, "commentCol"))
    if (length(commentCol) == 0) commentCol <- 0
    
    colClasses <- rep("NULL", nbCols)
    colClasses[dtCol] <- "character"
    colClasses[varCol] <- "numeric"
    if (commentCol > 0) colClasses[commentCol] <- "character"
    
    colNames <- paste("V", 1:nbCols)
    dtLab <- XML::xmlGetAttr(node, "dtLab")
    colNames[dtCol] <- dtLab
    colNames[varCol] <- varName
    if (commentCol > 0) colNames[commentCol] <- "comment"
    
    readData <- read.csv(file = path,
                         header = FALSE,
                         sep = XML::xmlGetAttr(node, "sep"),
                         skip = as.integer(XML::xmlGetAttr(node, "skip")),
                         colClasses = colClasses,
                         col.names= colNames)
    
    dtFormat <- XML::xmlGetAttr(node, "dtFormat")
    
    ## Now re-arrange datetime
    readData[ , dtCol] <- as.POSIXct(strptime(readData[ , dtCol],
                                              format = dtFormat, tz = TZ))
    if (commentCol == 0) {
        readData <- data.frame(readData,
                               comment = rep("", nrow(readData)),
                               stringsAsFactors = FALSE)
    }
    
    if (trace) {
        cat("events data read from file\n")
        if (nrow(readData) > 8L){
            print(head(readData, n = 4L))
            cat(" ... <more lines>\n")
            print(tail(readData, n = 4L))
        } else print(readData)
        
    }
    
    readData
    
}

## =============================================================================
##' Read a 'periods' file with description given in \code{node}
##'
##' @title Read a 'periods' file with description given in \code{node}
##'
##' @param node An XML node with attributes telling which columns are
##' \code{start}, \code{end}, and so on. The path of the file to read
##' (relative to \code{dir}) is given in the 'path' attribute of
##' 'node'.
##'
##' @param dir Path of the directory to read from.
##'
##' @param TZ Time zone.
##'
##' @param trace Level of verbosity.
##' 
##' @return A data frame with two \code{POSIXct} columns \code{start},
##' \code{end} a character column \code{comment} and possibly more.
##'
##' @author yves
##'
##' @note The \code{start} and \code{end} columns must be in the same
##' format given as \code{dtFormat} attribute of the node.  The file
##' read must be in the directory given in \code{dir}.
##' 
read.periodsFile <- function(node,
                             dir,
                             TZ = "GMT",
                             trace = 1) {
    
    stopifnot(requireNamespace("XML", quietly = TRUE)) 
    path <- file.path(dir, XML::xmlGetAttr(node, "path"))
    
    if (trace) cat("Reading the data file\n", path,"...\n")
    nbCols <- as.integer(XML::xmlGetAttr(node, "nbCols"))
    startCol <- as.integer(XML::xmlGetAttr(node, "startCol"))
    endCol <- as.integer(XML::xmlGetAttr(node, "endCol"))
    commentCol <- as.integer(XML::xmlGetAttr(node, "commentCol"))
    if (length(commentCol) == 0) commentCol <- 0
    
    dtFormat <- XML::xmlGetAttr(node, "dtFormat")
    
    colClasses <- rep("NULL", nbCols)
    colClasses[startCol] <- "character"
    colClasses[endCol] <- "character"
    if (commentCol > 0) colClasses[commentCol] <- "character"
    
    colNames <- paste("V", 1:nbCols)
    dtLAb <- XML::xmlGetAttr(node, "dtLab")
    colNames[startCol] <- "start"
    colNames[endCol] <- "end"
    if (commentCol > 0) colNames[commentCol] <- "comment"
    
    readData <- read.csv(file = path,
                         header = FALSE,
                         sep = XML::xmlGetAttr(node, "sep"),
                         skip = as.integer(XML::xmlGetAttr(node, "skip")),
                         colClasses = colClasses,
                         as.is = TRUE,
                         col.names = colNames)
    
    ## Now re-arrange datetimes
    readData[ , "start"] <- as.POSIXct(strptime(readData[ , "start"],
                                                format = dtFormat, tz = TZ))
    readData[ , "end"] <- as.POSIXct(strptime(readData[ , "end"],
                                              format = dtFormat, tz = TZ))
    
    if (trace) {
        cat("events data read from file\n")
        if (nrow(readData) > 8L){
            print(head(readData, n = 4L))
            cat(" ... <more lines>\n")
            print(tail(readData, n = 4L))
        } else print(readData)
        
    }
    
    readData
    
}

## =============================================================================
##' Read or Parse 'events' as specified in \code{node}
##'
##' For each child of \code{node} which is of type \code{events}
##' or \code{file}
##' 
##' @title Read or Parse 'events' as specified in \code{node}
##'
##' @param node An XML node.
##'
##' @param dir Path of the directory to read from.
##'
##' @param varName Variable name, as read in the ancestor of \code{node}.
##'
##' @param TZ Time zone.
##'
##' @param trace Level of verbosity.
##'
##' @return A data frame of events obtained by row binding of several
##' events data frames.
##'
##' @author yves
##'
##' @note The variable name must be given since it is an attribute of
##' an ancestor of \code{node}, and not of the node itself.
##' \code{node} should have children only of type \code{events} or
##' \code{file}. Anyway, other children will not be taken into
##' consideration. The node given in \code{node} formal is typically of
##' type "data".
##' 
readOrParse.events <- function(node,
                               dir,
                               varName = "Var",
                               TZ = "GMT",
                               trace = 1) {

    stopifnot(requireNamespace("XML", quietly = TRUE))
    
    events <- list()
    nevents <- 0
    ## file node???
    file.nodes   <- XML::xmlElementsByTagName(node, "file")
    
    if (length(file.nodes) > 1) {
        stop("at the time a 'data' can not contain more than one 'file' ",
             "child node")
    }
    
    if (length(file.nodes) > 0) {
        
        for (i in 1L:length(file.nodes) ) {
            nevents <- nevents + 1L
            events.i <- read.eventsFile(node = file.nodes[[i]],
                                        dir = dir,
                                        varName = varName,
                                        TZ = TZ,
                                        trace = trace)
            if (nevents > 1) events <- rbind(events, events.i)
            else  events <- events.i
        }
    }
    
    ## events node???
    events.nodes   <- XML::xmlElementsByTagName(node, "events")
    
    if (length(events.nodes) > 0) {
        for (i in 1L:length(events.nodes) ) {
            nevents <- nevents + 1
            events.i <- parse.eventsNode(node = events.nodes[[i]],
                                         varName = varName,
                                         TZ = TZ,
                                         trace = trace)
            if (nevents > 1) events <- rbind(events, events.i)
            else  events <- events.i
        }
    }
    
    ## returns a data.frame or NULL
    if (nevents) events
    else NULL
}

## =============================================================================
##' Read or Parse 'periods' as specified in an XML node
##' 
##' @title  Read or Parse 'periods' as specified in an XML node
##'
##' @param node An XML node.
##'
##' @param dir Path of the directory to read from.
##'
##' @param TZ Time zone.
##'
##' @param trace Level of verbosity.
##' 
##' @return A data frame of periods obtained by row binding of several
##' periods data frames.
##'
##' @author yves
##'
##' @note \code{node} should have children only of type \code{periods}
##' or \code{file}. Anyway, other children will not be taken into
##' consideration. The node given in \code{node} is typically of type
##' \code{data}.
##' 
readOrParse.periods <- function(node,
                                dir,
                                TZ = "GMT",
                                trace = 1) {
    
    stopifnot(requireNamespace("XML", quietly = TRUE))
    
    periods <- list()
    nperiods <- 0
    ## file node???
    file.nodes   <- XML::xmlElementsByTagName(node, "file")
    
    if (length(file.nodes) > 1)
        stop("at the time a 'data' can not contain more than one ",
             "'file' child node")
    
    if (length(file.nodes) > 0) {
        for (i in 1L:length(file.nodes) ) {
            nperiods <- nperiods + 1
            periods.i <- read.periodsFile(node = file.nodes[[i]],
                                          dir = dir,
                                          TZ = TZ,
                                          trace = trace)
            if (nperiods > 1) periods <- rbind(periods, periods.i)
            else  periods <- periods.i
        }
    }
    
    ## periods node???
    periods.nodes   <- XML::xmlElementsByTagName(node, "periods")
    
    if (length(periods.nodes) > 0) {
        for (i in 1L:length(periods.nodes) ) {
            nperiods <- nperiods + 1
            periods.i <- parse.periodsNode(node = periods.nodes[[i]],
                                           TZ = TZ,
                                           trace = trace)
            if (nperiods > 1) periods <- rbind(periods, periods.i)
            else  periods <- periods.i
        }
    }
    
    if (nperiods) periods
    else NULL
    
}
          
##=============================================================================
##' Read heterogenous data from an XML repository.    
##'
##' @title Read heterogenous data from an XML repository 
##'
##' @param name Name of the dataset to read.
##'
##' @param dir Path of the directory where files will be found.
##'
##' @param index Name of an index file which must be found in
##' \code{dir}. This file describes the structure of the datasets,
##' and one of them must have the name specified in \code{name}.
##' For some datasets, observations are read from csv files thet must
##' be located in \code{dir}.
##'
##' @param TZ Time zone.
##'
##' @param trace Level of verbosity.
##'
##' @return An object with class \code{"Rendata"}. This is a list with
##' \code{OTinfo} \code{OTdata}. It also can have elements
##' \code{OTSinfo} and \code{OTSdata} if \code{OTS} blocks are given
##' in the file or \code{MAXinfo} and \code{MAXdata} if \code{MAX}
##' blocks are given in the file.
##'
##' @author yves
##' 
readXML <- function(name,
                    dir,
                    index = "index.xml",
                    TZ = "GMT",
                    trace = 0) {
    
    stopifnot(requireNamespace("XML", quietly = TRUE))
    
    ##-------------------------------------------------------------------------
    ## check that the file exists and can be read
    ##-------------------------------------------------------------------------

    if ((file.access(names = dir, mode = 0)!=0) || (!file.info(dir)$isdir)) {
        stop("dir =", dir, "is not an existing directory")
    }
    file <- file.path(dir, index)
    
    if (trace) cat("Reading 'index' file\n", file, "...\n")
    
    if (file.access(names = file, mode = 4) != 0) {
        if (file.access(names = file, mode = 0) != 0) {
            stop("'index' file not found")
        } else stop("not allowed to open file")
    }
    
    arbre   <- XML::xmlTreeParse(file)
    racine  <- XML::xmlRoot(arbre)
    
    dataset.node   <- XML::xmlElementsByTagName(racine, "dataset")
    dataset.names <- sapply(dataset.node, XML::xmlGetAttr, "name")
    
    if (trace) {
        cat("datasets declared\n")
        print(dataset.names)
    }
    
    m <- match(dataset.names, name)
    
    if (trace) {
        cat("matching...\n")
        print(m)
    }
    
    nm <- sum(m == 1, na.rm = TRUE)
    
    if (nm==0) stop("no correct dataset")
    else if (nm>1) stop("several datasets match 'name'")
    
    ind <- (1:length(dataset.node))[!is.na(m) & m==1]
    dataset.node <- dataset.node[[ind]]
    
    ##=========================================================================
    ## Information about the dataset
    ##=========================================================================
    
    info <- list()
    for (attr  in c("name", "shortLab", "longLab", "varName",
                    "varShortLab", "varUnit")) {
        info[[attr]] <- XML::xmlGetAttr(dataset.node, attr)
    }
    
    describe.node   <- XML::xmlElementsByTagName(dataset.node, "describe")
    if (length(describe.node) == 0) describe.node <- NULL
    else describe.node   <- as(describe.node[[1]], "character")
    
    ##========================================================================
    ## "Over Threshold" data. The dataset must contain EXACTLY ONE
    ## OTdata node.
    ##========================================================================
    
    OTdata.nodes   <- XML::xmlElementsByTagName(dataset.node, "OTdata")
    if ( (length(OTdata.nodes) > 1) || (length(OTdata.nodes) == 0 ) )
        stop("a 'dataset' node must contain exactly one 'OTdata' child node")
    
    OTinfo <- list()
    for (attr in c("start", "end")) 
        OTinfo[[attr]] <- as.POSIXct(XML::xmlGetAttr(OTdata.nodes[[1]], attr),
                                     tz = TZ)
    
    for (attr in c("effDuration", "threshold"))
        OTinfo[[attr]] <- as.numeric(XML::xmlGetAttr(OTdata.nodes[[1]], attr))
    
    data.nodes   <- XML::xmlElementsByTagName(OTdata.nodes[[1]], "data")
    if ( length(data.nodes) != 1 ) 
        stop("an 'OTdata' node must contain exactly one 'data' child node")
    
    OTdata <- readOrParse.events(node = data.nodes[[1]],
                                 dir = dir, varName = info$varName,
                                 TZ = TZ,
                                 trace = trace)
    
    ##print(OTdata)
    
    if (any(OTdata$date < OTinfo$start)) {
        warning("'OTdata' contain events with date < OTinfo$start.",
                " They are removed.")
        OTdata <- OTdata[OTdata$date >= OTinfo$start, ]
    }
    
    if (any(OTdata[ , "date"] > OTinfo$end)) {
        warning("'OTdata' contain events with date > OTinfo$end.",
                " They are removed. ")
        OTdata <- OTdata[OTdata$date <= OTinfo$end, ]
    }
    
    if (any(OTdata[ , info$varName] < OTinfo$threshold)) {
        warning("'OTdata' column '", info$varName,
                "' contains a value < info$threshold")
    }
    
    ##=========================================================================
    ## Missing periods for OTdata We must find 'start' and 'end' for
    ## the successive periods, and caution is needed for successive
    ## periods possibly collapsing.
    ##=========================================================================
    
    missing.nodes   <- XML::xmlElementsByTagName(OTdata.nodes[[1]], "missing")
    
    if( length(missing.nodes) > 1)  
        stop("an \"OTdata\" node can not have more than one child  \"missing\"")
    
    if( length(missing.nodes) == 1) {
        
        OTmissing <- readOrParse.periods(node = missing.nodes[[1]],
                                         dir = dir,
                                         TZ = TZ,
                                         trace = trace)
        
        OTmissing <- cleanInt(data = OTmissing,
                              start = OTinfo$start,
                              end = OTinfo$end,
                              TZ = TZ)
        
        ## warn if some OT events fall into an OTmissing period
        inGaps <- inInts(x = OTdata$date,
                         periods = OTmissing)
        
        if (any(inGaps)) {
            warning(sum(inGaps),"'OTdata' events are in 'OTmissing' periods.",
                    " They are removed.")
            OTdata <- OTdata[inGaps, ]
        }
        
        ## Check effective duration
        noskip <- skip2noskip(skip = OTmissing,
                              start = OTinfo$start,
                              end = OTinfo$end)
        
        checkDuration <-
            sum(as.numeric(difftime(noskip$end, noskip$start, units = "day") /
                               365.25))
        
        if (abs(checkDuration - OTinfo$effDuration) > 0.1)
            warning("'effDuration' attributes does not agrre with computation")
        
        if (trace) cat("checkDuration", checkDuration, "effDuration",
                       OTinfo$effDuration, "\n")
        
    } else {
        checkDuration <-
            sum(as.numeric(difftime(OTinfo$end, OTinfo$start, units = "day") /
                               365.25))
        
        if (abs(checkDuration - OTinfo$effDuration) > 0.1)
            warning("'effDuration' attributes does not agree with computation")
        
        if (trace) cat("checkDuration", checkDuration, "effDuration",
                       OTinfo$effDuration, "\n")
        
        OTmissing <- NULL
    }
    
    ##=========================================================================
    ## Prepare results list
    ##=========================================================================
    
    MyList <- list(info = info,
                   describe = describe.node,
                   OTinfo = OTinfo,
                   OTdata = OTdata,
                   OTmissing = OTmissing)
    
    ##=========================================================================
    ## MAX data nodes
    ##=========================================================================
    
    MAXdata.nodes   <- XML::xmlElementsByTagName(dataset.node, "MAXdata")
    MAXinfo <- list()
    
    if (length(MAXdata.nodes)) {
        
        if (trace) {
            cat("Processing", length(MAXdata.nodes),
                "'MAX' (historical) data\n")
        }
        
        MAXdata <- list()
        
        for (i in 1:length(MAXdata.nodes)){
            
            node <-  MAXdata.nodes[[i]]
            
            if (i == 1)  MAXinfo[["block"]] <- 1
            else MAXinfo[["block"]] <- c(MAXinfo[["block"]], i)
            
            for (attr  in c("start", "end")) {
                if (i == 1)  MAXinfo[[attr]] <- XML::xmlGetAttr(node, attr)
                else {
                    MAXinfo[[attr]] <- c(MAXinfo[[attr]],
                                         XML::xmlGetAttr(node, attr))
                }
            }
            
            ##-----------------------------------------------------------------
            ## data (events) Note that there can be no events in some
            ## case.  Then there is a 'data' node but with no 'events'
            ## nor 'file' child.
            ## ----------------------------------------------------------------
            
            data.nodes   <- XML::xmlElementsByTagName(node, "data")
            
            if(length(data.nodes) != 1L) {  
                stop("a \"MAXdata\" node must have exactly one child \"data\"")
            }
            MAXdata.i <- readOrParse.events(data.nodes[[1]],
                                            dir = dir, varName = info$varName,
                                            TZ = TZ,
                                            trace = trace)
            
            ri <- 0
            if (!is.null(MAXdata.i)) {
                ri <- nrow(MAXdata.i)
                MAXdata.i <- data.frame(block = rep(i, times = ri),
                                        MAXdata.i,
                                        stringsAsFactors = FALSE)
            } 
            
            if (i == 1) MAXdata <- MAXdata.i
            else if (ri) MAXdata <- rbind(MAXdata, MAXdata.i)
            
            if (i == 1)  MAXinfo[["r"]] <- ri
            else MAXinfo[["r"]] <- c(MAXinfo[["r"]], ri)
            
            if (is.null(MAXdata.i)) {
                warning("'MAXdata' ", i, " contains no events!")
            }
            
            dc <- MAXdata.i$date[!is.na(MAXdata.i$date)]
            
            if (!is.null(MAXdata.i) &&
                any(dc < as.POSIXct(MAXinfo$start[i], tz = TZ))) {
                warning("'MAXdata' ", i, " contains an event before",
                        " 'start' given in 'MAXinfo'!")
            }
            if (!is.null(MAXdata.i) &&
                any(dc > as.POSIXct(MAXinfo$end[i], tz = TZ))) {
                warning("'MAXdata' ", i, " contains an event after",
                        " 'end' given in 'MAXinfo'!")
            }
            
        }
        
        ## MAXinfo <- as.data.frame(MAXinfo)
        MAXinfo$start <- as.POSIXct(MAXinfo$start, tz = TZ)
        MAXinfo$end <- as.POSIXct(MAXinfo$end, tz = TZ)
        MAXinfo$r <- as.integer(MAXinfo$r)
        ## Only if missing
        dur <- difftime(MAXinfo$end, MAXinfo$start, units= "day") / 365
        dur <- round(as.numeric(dur),  digits = 2)
        MAXinfo <- data.frame(start = MAXinfo$start,
                              end =  MAXinfo$end,
                              duration = dur)
        
        MyList$MAXinfo <- MAXinfo
        MyList$MAXdata <- MAXdata
        
    } else {
        if (trace) cat("No MAX (historical) data\n")
    }
    
    ##========================================================================
    ## OTS data nodes
    ##========================================================================
    
    OTSdata.nodes   <- XML::xmlElementsByTagName(dataset.node, "OTSdata")
    
    if (length(OTSdata.nodes)) {
        
        if (trace) cat("Processing", length(OTSdata.nodes),
                       "'OTS' (historical) data\n")
        
        OTSinfo <- list()
        
        for (i in 1L:length(OTSdata.nodes)) {
            
            node <-  OTSdata.nodes[[i]]

            if (i == 1)  {
                OTSinfo[["block"]] <- 1
            }
            
            if (i == 1)  OTSinfo[["block"]] <- 1
            else OTSinfo[["block"]] <- c(OTSinfo[["block"]], i)
            
            for (attr  in c("start", "end", "threshold")) {
                if (i == 1)  OTSinfo[[attr]] <- XML::xmlGetAttr(node, attr)
                else {
                    OTSinfo[[attr]] <- c(OTSinfo[[attr]],
                                         XML::xmlGetAttr(node, attr))
                }
            }
            
            ##-----------------------------------------------------------------
            ## data (events) Note that there can be no events in some
            ## case.  Then there is a 'data' node but with no 'events'
            ## nor 'file' child.
            ## ----------------------------------------------------------------
            data.nodes   <- XML::xmlElementsByTagName(node, "data")
            
            if( length(data.nodes) != 1)  
                stop("a \"MAXdata\" node must have exactly one child \"data\"")
            
            OTSdata.i <- readOrParse.events(data.nodes[[1]],
                                            dir = dir, varName = info$varName,
                                            TZ = TZ,
                                            trace = trace)
            ri <- 0
            if (!is.null(OTSdata.i)) {
                ri <- nrow(OTSdata.i)
                OTSdata.i <- data.frame(block = rep(i, times = ri),
                                        OTSdata.i,
                                        stringsAsFactors = FALSE)
            } 
            
            if (i == 1) OTSdata <- OTSdata.i
            else if (ri) OTSdata <- rbind(OTSdata, OTSdata.i)
            
            if (i == 1)  OTSinfo[["r"]] <- ri
            else OTSinfo[["r"]] <- c(OTSinfo[["r"]], ri)

            if (!is.null(OTSdata.i)) {
                dc <- OTSdata.i$date[!is.na(OTSdata.i$date)]

                if (any(dc < as.POSIXct(OTSinfo$start[i]))) {
                    warning("'OTSdata' ", i, " contains an event ",
                            "before 'start' given in 'OTSinfo'")
                }
                if (any(dc > as.POSIXct(OTSinfo$end[i]))) {
                    warning("'OTSdata' ", i, " contains an event ",
                            "after 'end' given in 'OTSinfo'")
                }
            }
            
        }
        
        ## OTSinfo$start <- as.POSIXct(OTSinfo$start)
        ## OTSinfo$end <- as.POSIXct(OTSinfo$end)
        dur <- difftime(OTSinfo$end, OTSinfo$start, units= "day") / 365
        dur <- round(as.numeric(dur),  digits = 2)
        
        OTSinfo <- data.frame(start = as.POSIXct(OTSinfo$start, tz = TZ),
                              end = as.POSIXct(OTSinfo$end, tz = TZ),
                              duration = dur,
                              threshold = as.numeric(OTSinfo$threshold),
                              r = as.integer(OTSinfo$r))
        
        if (any(OTSinfo$threshold < OTinfo$threshold) )
            warning("'OTSdata' specified with a threshold < OTinfo$threshold")
        
        MyList$OTSinfo <- OTSinfo
        MyList$OTSdata <- OTSdata
        
    } else {
        if (trace) cat("No OTS (historical) data\n")
    }
    
    
    attr(MyList, "class") <- "Rendata"
    MyList
    
}
