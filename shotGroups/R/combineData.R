combineData <-
function(DFs) {
    if(!is.list(DFs)) { stop("DFs must be a list") }
    if(!all(vapply(DFs, is.data.frame, logical(1)))) {
        stop("DFs must be a list of data frames")
    }

    ## shared set of variable names
    varNameL <- lapply(DFs, names)       # list of variable names
    varNames <- tolower(Reduce(intersect, varNameL))

    ## check if data frames contain required variables
    wantsGrp <- "group"                  # useful
    wantsDst <- "distance"               # useful
    hasGrp   <- wantsGrp %in% varNames   # useful ones we have
    hasDst   <- wantsDst %in% varNames   # useful ones we have

    if(!all(hasGrp)) {
        warning(c("At least one data frame is missing variable\n",
                  paste(wantsGrp[!hasGrp], collapse=" "),
                  "\nGroup is set to 1"))
    }

    ## if variable group does not exist, is NA or empty - set
    setGroup <- function(x) {
        if(!("group" %in% names(x))) {
            x$group <- 1
        } else if(all(x$group == "") || all(is.na(x$group))) {
            x$group <- 1
        }

        x
    }

    DFs <- lapply(DFs, setGroup)

    if(!all(hasDst)) {
        warning(c("At least one file is missing variable\n",
                  paste(wantsDst[!hasDst], collapse=" "),
                  "\nthat may be required later by analysis functions"))
    }

    ## make sure each data frame has either X, Y or Point.X, Point.Y
    replaceXY <- function(x) {
        x        <- setNames(x, tolower(names(x))) # convert to lower case
        dfNames  <- names(x)
        needsXY1 <- c("point.x", "point.y")  # coordinates must have this name
        needsXY2 <- c("x", "y")              # or this
        hasXY1   <- needsXY1 %in% dfNames
        hasXY2   <- needsXY2 %in% dfNames

        if(!xor(all(hasXY1), all(hasXY2))) { # not (either X, Y or Point.X, Point.Y)
            stop("Coordinates must be named either X, Y or Point.X, Point.Y")
        }

        if(("z" %in% dfNames) && ("point.z" %in% dfNames)) {
            stop("Coordinates must be named either Z or Point.Z")
        }

        ## if X, Y -> rename to Point.X, Point.Y
        if(all(hasXY2)) {
            dfNames[dfNames %in% "x"] <- "point.x"
            dfNames[dfNames %in% "y"] <- "point.y"
            dfNames[dfNames %in% "z"] <- "point.z"
            warning("Variables X, Y were renamed to Point.X, Point.Y")
            names(x) <- dfNames
        }

        x
    }

    DFs <- lapply(DFs, replaceXY)

    ## add new group variable that is more descriptive
    setGroupVerbose <- function(x) {
        ## if project title is available -> use it
        ## if not -> use file name
        ## if ammunition is available -> use it
        groupVA <- if(!is.null(x$project.title)) {
            if(!all(is.na(x$project.title)) && !all(x$project.title == "")) {
                x$project.title
            } else {
                x$file
            }
        } else {
            x$file
        }

        groupVB <- if(!is.null(x$ammunition)) {
            if(!all(is.na(x$ammunition)) && !all(x$ammunition == "")) {
                x$ammunition
            } else {
                ""
            }
        } else {
            ""
        }

        #groupVerb <- paste(x$group, groupVA, groupVB, sep="_")
        groupVerb <- paste(groupVA, groupVB, sep="_")

        ## trim leading/trailing _, collapse __ to _, replace " " with _
        groupVerb   <-  sub("_$", "",  groupVerb)
        groupVerb   <-  sub("^_", "",  groupVerb)
        groupVerb   <- gsub("__", "_", groupVerb)
        groupVerb   <- gsub("[[:blank:]]+", "_", groupVerb)
        x$groupVerb <- groupVerb
        x
    }

    DFs <- lapply(DFs, setGroupVerbose)

    ## restrict data frames to shared variables variables
    varsNow <- Reduce(intersect, lapply(DFs, names))  # shared set of variables
    DFrestr <- lapply(DFs, function(x) x[, varsNow])  # select only these
    nObs    <- vapply(DFrestr, nrow, integer(1))      # number of observations in each data frame
    DFall   <- do.call("rbind", DFrestr)     # combine data frames
    rownames(DFall) <- NULL                  # remove row names

    ## add new factor Origin for coding original file
    origin <- factor(rep(seq_along(DFs), nObs))

    ## add new factor series for coding Groups as a consecutive number over files
    ## first a factor with alphabetically ordered levels
    orgSer <- droplevels(interaction(origin, DFall$group))

    ## convert orgSer to a factor with consecutively numbered levels + description
    runs         <- rle(as.character(orgSer))
    runs$values  <- seq_along(runs$values)
    seriesNum    <- inverse.rle(runs)
    DFall$series <- factor(paste(seriesNum, DFall$groupVerb, sep="_"),
                           labels=unique(paste(seriesNum, DFall$groupVerb, sep="_")))

    ## convert orgSer to a factor with consecutively numbered levels
    DFall$seriesNum <- seriesNum

    return(DFall)
}
