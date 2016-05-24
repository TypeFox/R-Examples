# Automatically generated from the noweb directory
tmerge <- function(data1, data2, id, ..., tstart, tstop, options) {
    Call <- match.call()
    # The function wants to recognize special keywords in the
    #  arguments, so define a set of functions which will be used to
    #  mark objects
    new <- new.env(parent=parent.frame())
    assign("tdc", function(time, value=NULL) {
        x <- list(time=time, value=value); 
        class(x) <- "tdc"; x},
           envir=new)
    assign("cumtdc", function(time, value=NULL) {
        x <- list(time=time, value=value); 
        class(x) <-"cumtdc"; x},
           envir=new)
    assign("event", function(time, value=NULL, censor=NULL) {
        x <- list(time=time, value=value, censor=censor); 
        class(x) <-"event"; x},
           envir=new)
    assign("cumevent", function(time, value=NULL, censor=NULL) {
        x <- list(time=time, value=value, censor=censor); 
        class(x) <-"cumevent"; x},
           envir=new)

    if (missing(data1) || missing(data2) || missing(id)) 
        stop("the data1, data2, and id arguments are required")
    if (!inherits(data1, "data.frame")) stop("data1 must be a data frame")
    
    tmerge.control <- function(idname="id", tstartname="tstart", tstopname="tstop",  
                               delay =0, na.rm=TRUE, ...) {
        extras <- list(...)
        if (length(extras) > 0) 
            stop("unrecognized option(s):", paste(names(extras), collapse=', '))
        if (length(idname) != 1 || make.names(idname) != idname)
            stop("idname option must be a valid variable name")
        if (!is.null(tstartname) && 
            (length(tstartname) !=1 || make.names(tstartname) != tstartname))
            stop("tstart option must be NULL or a valid variable name")
        if (length(tstopname) != 1 || make.names(tstopname) != tstopname)
            stop("tstop option must be a valid variable name") 
        if (length(delay) !=1 || !is.numeric(delay) || delay < 0)
            stop("delay option must be a number >= 0")
        if (length(na.rm) !=1 || ! is.logical(na.rm))
            stop("na.rm option must be TRUE or FALSE")
       list(idname=idname, tstartname=tstartname, tstopname=tstopname, delay=delay, na.rm=na.rm)
    }

    tname <- attr(data1, "tname")
    firstcall <- is.null(tname)  #first call to the function
    if (!firstcall && any(is.null(match(unlist(tname), names(data1)))))
        stop("data1 does not match its own tname attribute")

    if (!missing(options)) {
        if (!is.list(options)) stop("options must be a list")
        if (!is.null(tname)) {
            # If an option name matches one already in tname, don't confuse
            #  the tmerge.control routine with duplicate arguments
            temp <- match(names(options), names(tname), nomatch=0)
            topt <- do.call(tmerge.control, c(options, tname[temp==0]))
            if (any(temp >0)) {
                # A variable name is changing midstream, update the
                # variable names in data1
                varname <- tname[c("idname", "tstartname", "tstopname")]
                temp2 <- match(varname, names(data1))
                names(data1)[temp2] <- varname
            }
        }
        else topt <- do.call(tmerge.control, options)
    }
    else if (length(tname)) topt <- do.call(tmerge.control, tname)
    else  topt <- tmerge.control()

    # id, tstart, tstop are found in data2
    if (missing(id)) stop("the id argument is required")
    if (missing(data1) || missing(data2))
        stop("two data sets are required")
    id <- eval(Call[["id"]], data2, enclos=emptyenv()) #don't find it elsewhere
    if (is.null(id)) stop("the id variable is null")

    if (firstcall) {
        if (!missing(tstop)) {
             tstop <-  eval(Call[["tstop"]],  data2)
             if (length(tstop) != length(id))
                 stop("tstop and id must be the same length")
             # The neardate routine will check for legal tstop data type
          }
        if (!missing(tstart)) {
            tstart <- eval(Call[["tstart"]], data2)
            if (length(tstart)==1) tstart <- rep(tstart, length(id))
            if  (length(tstart) != length(id))        
                stop("tstart and id must be the same length")
            if (any(tstart >= tstop))
                stop("tstart must be < tstop")
             }
    }
    else {
        if (!missing(tstart) || !missing(tstop))
            stop("tstart and tstop arguments only apply to the first call")
    }
    # grab the... arguments
    notdot <- c("data1", "data2", "id", "tstart", "tstop", "options")
    dotarg <- Call[is.na(match(names(Call), notdot))]
    dotarg[[1]] <- as.name("list")  # The as-yet dotarg arguments
    if (missing(data2)) args <- eval(dotarg, envir=new)
    else  args <- eval(dotarg, data2, enclos=new)
        
    argclass <- sapply(args, function(x) (class(x))[1])
    argname <- names(args)
    if (any(argname== "")) stop("all additional argments must have a name")
           
    check <- match(argclass, c("tdc", "cumtdc", "event", "cumevent"))
    if (any(is.na(check)))
        stop(paste("argument(s)", argname[is.na(check)], 
                       "not a recognized type"))
    # The tcount matrix is useful for debugging
    tcount <- matrix(0L, length(argname), 8)
    dimnames(tcount) <- list(argname, c("early","late", "gap", "within", 
                                        "boundary", "leading", "trailing",
                                        "tied"))
    tevent <- attr(data1, "tevent") # event type variables
    tcens  <- attr(data1, "tcensor")# censor code for variables
    if (is.null(tcens)) tcens <- vector('list', 0)
    newdata <- data1 #make a copy
    if (firstcall) {
        # The line below finds id, tstop, and tstart variables in data1
        indx <- match(c(topt$idname, topt$tstartname, topt$tstopname), names(data1), 
                      nomatch=0)
        if (any(indx[2:3]>0) && FALSE) {  # warning currently turned off. Be chatty?
            overwrite <- c(topt$tstartname, topt$tstopname)[indx[2:3]]
            warning("overwriting data1 variables", paste(overwrite, collapse=' '))
            }
        
        if (indx[1] == 0) {
            # the topt$id variable name is not in data1.  Deal with the
            # fairly common case that data1 == data2.  If data1 has the
            # variable used as 'id' in data2, add it to data1 under the
            # new name
            temp <- as.character(Call[["id"]])
            if (is.name(Call[["id"]]) && !is.na(match(temp, names(data1)))) {
                data1[[topt$idname]] <- data1[[temp]]
                baseid <- data1[[temp]]
                }
            else stop("id variable not found in data1")
        }
        else baseid <- data1[[indx[1]]]

        if (any(duplicated(baseid))) 
            stop("for the first call (that establishes the time range) data1 must have no duplicate identifiers")

        if (length(baseid)== length(id) && all(baseid == id)) newdata <- data1
        else {
            indx2 <- match(id, baseid)
            if (any(is.na(indx2)))
                stop("'id' has values not in data1")
            newdata <- data1[indx2,]
            }
        if (missing(tstop)) { # case 2
            if (length(argclass)==0 || argclass[1] != "event")
                stop("neither a tstop argument nor an initial event argument was found")
            tstop <- args[[1]][[1]]
            }
         
        # at this point newdata and data2 are in the same order, same # rows
        if (any(is.na(tstop))) 
            stop("missing time value, when that variable defines the span")
        if (missing(tstart)) tstart <- rep(0, length(id))
        if (any(tstart >= tstop)) 
            stop("tstart must be > tstop")
        newdata[[topt$tstartname]] <- tstart
        newdata[[topt$tstopname]] <- tstop
        if (any(duplicated(id))) {
            # sort by time within id
            indx1 <- match(id, unique(id))
            newdata <- newdata[order(indx1, tstop),]
        }
        n <- nrow(newdata)
        temp <- newdata[[topt$idname]]
        if (any(tstart >= tstop)) stop("tstart must be < tstop")
        if (any(newdata$tstart[-n] > newdata$tstop[-1] &
                temp[-n] == temp[-1]))
            stop("there are overlapping time intervals")
    }
    else { #not a first call
        if (any(is.na(match(id, data1[[topt$idname]]))))
            stop("id values were found in data2 which are not in data1")
    }
    saveid <- id
    for (ii in seq(along.with=args)) {
        argi <- args[[ii]]
        baseid <- newdata[[topt$idname]]
        dstart <- newdata[[topt$tstartname]]
        dstop  <- newdata[[topt$tstopname]]
        argcen <- argi$censor
        
        # if an event time is missing then skip that obs
        etime <- argi$time
        if (length(etime) != length(saveid))
            stop("argument ", argname[ii], " is not the same length as id")
       if (!is.null(argi$value)) {
           if (length(argi$value) != length(saveid))
                stop("argument", argname[ii], "is not the same length as id")
            if (topt$na.rm) keep <- !(is.na(etime) | is.na(argi$value))
            else keep <- !is.na(etime)
            if (!all(keep)) {
                etime <- etime[keep]
                argi$value <- argi$value[keep]
                }
            }
        else {
          keep <- !is.na(etime)
          etime <- etime[keep]
          }
        id <- saveid[keep]
        
        # For an event or cumevent, one of the later steps becomes much
        #  easier if we sort the new data by id and time
        indx <- order(id, etime)
        id <- id[indx]
        etime <- etime[indx]
        if (!is.null(argi$value))
            yinc <- argi$value[indx]
            
        # indx1 points to the closest start time in the baseline data (data1)
        #  that is <= etime.  indx2 to the closest end time that is >=etime.
        # If etime falls into a (tstart, tstop) interval, indx1 and indx2
        #   will match
        # If the "delay" argument is set and this event is of type tdc, then
        #   move any etime that is after the entry time for a subject.
        if (topt$delay >0 && argclass[ii] %in% c("tdc", "cumtdc")) {
            mintime <- tapply(dstart, baseid, min)
            index <- match(id, names(mintime))
            etime <- ifelse(etime <= mintime[index], etime, etime+ topt$delay)
        }
        
        indx1 <- neardate(id, baseid, etime, dstart, best="prior")
        indx2 <- neardate(id, baseid, etime, dstop, best="after")

        # The event times fall into one of 5 categories
        #   1. Before the first interval
        #   2. After the last interval
        #   3. Outside any interval but with time span, i.e, it falls into
        #       a gap in follow-up
        #   4. Strictly inside an interval (does't touch either end)
        #   5. Inside an interval, but touching.
        itype <- ifelse(is.na(indx1), 1,
                        ifelse(is.na(indx2), 2, 
                               ifelse(indx2 > indx1, 3,
                                      ifelse(etime== dstart[indx1] | 
                                             etime== dstop[indx2], 5, 4))))

        # Subdivide the events that touch on a boundary
        #  1: intervals of (a,b] (b,d], new count at b  "tied edge"
        #  2: intervals of (a,b] (c,d] with c>b, new count at c, "front edge"
        #  3: intervals of (a,b] (c,d] with c>b, new count at b, "back edge"
        #
        subtype <- ifelse(itype!=5, 0, 
                          ifelse(indx1 == indx2+1, 1,
                                 ifelse(etime==dstart[indx1], 2, 3)))
        tcount[ii,1:7] <- table(factor(itype+subtype, levels=c(1:4, 6:8)))

        # count ties.  id and etime are not necessarily sorted
        tcount[ii,8] <- sum(tapply(etime, id, function(x) sum(duplicated(x))))
        # Look to see if this term has one or two arguments.  If one arg
        #  then the increment is 1, else it is the second arg.
        if (!is.null(argi$value)) yinc <- argi$value
        else yinc <- rep(1.0, length(etime))
        indx4 <- which(itype==4)
        n4 <- length(indx4)
        if (n4 > 0) {
            icount <- tapply(etime[indx4], indx1[indx4], function(x) sort(unique(x)))
            n.add <- sapply(icount, length)  #number of rows to add
            
            # expand the data 
            irep <- rep.int(1L, nrow(newdata))
            erow <- unique(indx1[indx4])   # which rows in newdata to be expanded
            irep[erow] <- 1+ n.add # number of rows in new data
            jrep <- rep(1:nrow(newdata), irep)  #stutter the duplicated rows
            newdata <- newdata[jrep,]  #expand it out
            dstart <- dstart[jrep]
            dstop <-  dstop[jrep]

            #fix up times
            nfix <- length(erow)
            temp <- vector("list", nfix)
            iend <- (cumsum(irep))[irep >1]  #end row of each duplication set
            for (j in 1:nfix) temp[[j]] <-  -(seq(n.add[j] -1, 0)) + iend[j]
            newrows <- unlist(temp)
            dstart[newrows] <- dstop[newrows-1] <- unlist(icount)
            newdata[[topt$tstartname]] <- dstart
            newdata[[topt$tstopname]]  <- dstop
            for (ename in tevent) newdata[newrows-1, ename] <- tcens[[ename]]

            # refresh indices
            baseid <- newdata[[topt$idname]]
            indx1 <- neardate(id, baseid, etime, dstart, best="prior")
            indx2 <- neardate(id, baseid, etime, dstop, best="after")
            subtype[itype==4] <- 1  #all the "insides" are now on a tied edge
            itype[itype==4]   <- 5  
        }
        # add it in
        if (argclass[ii] %in% c("cumtdc", "cumevent")) {
            if (!is.numeric(yinc)) stop("invalid increment for cumtdc or cumevent")
            yinc <- unlist(tapply(yinc, id, cumsum))
            }

        newvar <- newdata[[argname[ii]]]  #does the variable exist? 
        if (argclass[ii] %in% c("event", "cumevent")) {
            if (is.null(newvar)) {
                if (is.factor(yinc)) 
                    newvar <- factor(rep(levels(yinc)[1], nrow(newdata)),
                                     levels(yinc))
                else if (is.numeric(yinc)) newvar <- rep(0L, nrow(newdata))
                else stop("invalid value for a status variable")
            }
            keep <- (subtype==1 | subtype==3) # all other events are thrown away
            newvar[indx2[keep]] <- yinc[keep]
            
            if (!(argname[ii] %in% tevent)) {
                tevent <- c(tevent, argname[[ii]])
                if (is.factor(yinc)) tcens <- c(tcens, levels(yinc)[1])
                else tcens <- c(tcens, 0)
                names(tcens) <- tevent
                }
        }
        else {
            keep <- itype != 2  # changes after the last interval are ignored
            indx <- ifelse(subtype==1, indx1, 
                           ifelse(subtype==3, indx2+1L, indx2))

            if (is.null(newvar)) {
                if (is.null(argi$value)) newvar <- rep(0.0, nrow(newdata))
                else newvar <- rep(NA_real_, nrow(newdata))
                }
            # id can be any data type; feed integers to the C routine
            storage.mode(yinc) <- storage.mode(dstop) <- "double"
            storage.mode(newvar) <- storage.mode(etime) <- "double"
            newvar <- .Call("tmerge", match(baseid, baseid), dstop, newvar, 
                            match(id, baseid)[keep], etime[keep], 
                            yinc[keep], indx[keep])
        }

        newdata[[argname[ii]]] <- newvar
    }
    attr(newdata, "tname") <- topt[c("idname", "tstartname", "tstopname")]
    attr(newdata, "tcount") <- rbind(attr(data1, "tcount"), tcount)
    if (length(tevent)) {
        attr(newdata, "tevent") <- tevent
        attr(newdata, "tcensor" ) <- tcens
        }
    row.names(newdata) <- NULL
    class(newdata) <- c("data.frame")
    newdata
}
