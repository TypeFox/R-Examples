## read files that are comma or whitespace-delimited
## variable names must not have spaces
readDataMisc <-
function(fPath=".", fNames, fPat, combine=TRUE) {
    files <- getFileNames(fPath, fNames, fPat)
    if(length(files) < 1L) { stop("No files were selected") }

    ## determine whether file is csv or whitespace delimited
    isCSV <- function(ext, nFields) {
        csv <- if(tolower(ext) == "csv") {        # explicit file extension
            TRUE
        } else if((length(unique(nFields)) == 1L) && (nFields[1] >= 2)) {
            ## same number of comma-separated fields in all rows
            ## and at least two comma-separated fields
            TRUE
        } else {
            FALSE
        }
    }

    ## read a single file, possibly csv, possibly whitespace delimited
    readMe <- function(f) {
        pieces  <- strsplit(f, "\\.")[[1]]
        ext     <- tolower(pieces[length(pieces)]) # file extensions
        nFields <- count.fields(f, sep=",")

        ## choose appropriate function to read in file
        readFun <- if(isCSV(ext, nFields)) {
            read.csv
        } else {
            read.table
        }

        DF      <- readFun(f, header=TRUE)
        DF$file <- basename(tools::file_path_sans_ext(f))  # add filename
        setNames(DF, tolower(names(DF))) # convert var names to lower case
    }

    ## read multiple files, some possibly csv, some possibly whitespace delimited
    DFs <- lapply(files, readMe)
    names(DFs) <- paste("file", seq_along(DFs), sep="")  # name them

    ## build shared set of variable names
    varNameL <- lapply(DFs, names)           # list of variable names
    varNames <- Reduce(intersect, varNameL)  # intersection of all var names

    ## make sure that the data frames all have the correct variables
    ## remove dots from variable names
    wants   <- c("distance", "group", "aimx", "aimy")  # useful
    vnNoDot <- vapply(strsplit(varNames, "\\."), function(y) { paste0(y, collapse="") }, character(1))
    has     <- wants %in% vnNoDot

    if(!all(has)) {
        warning(c("At least one file is missing variable(s)\n",
            paste(wants[!has], collapse=" "),
            "\nthat may be required later by analysis functions"))
    }

    ## make sure each data frame has either X, Y or Point.X, Point.Y
    replaceXY <- function(x) {
        dfNames  <- names(x)
        needsXY1 <- c("point.x", "point.y")  # coordinates must have this name
        needsXY2 <- c("x", "y")              # or this
        needsXY3 <- c("shotx", "shoty")      # or this (Taran)
        hasXY1   <- needsXY1 %in% dfNames
        hasXY2   <- needsXY2 %in% dfNames
        hasXY3   <- needsXY3 %in% dfNames

        ## exactly one of X, Y - Point.X, Point.Y - ShotX, ShotY
        if(sum(c(all(hasXY1), all(hasXY2), all(hasXY3))) != 1) {
            stop("Coordinates must be named X, Y - Point.X, Point.Y - or ShotX, ShotY")
        }

        if(sum(c("z" %in% dfNames, "point.z" %in% dfNames, "shotz" %in% dfNames)) > 1L) {
            stop("Coordinates must be named Z, Point.Z, or ShotZ")
        }

        ## if X, Y, Z or ShotX, ShotY, ShotZ -> rename to Point.X, Point.Y, Point.Z
        if(all(hasXY2) | all(hasXY3)) {
            dfNames[dfNames %in% c("x", "shotx")] <- "point.x"
            dfNames[dfNames %in% c("y", "shoty")] <- "point.y"
            dfNames[dfNames %in% c("z", "shotz")] <- "point.z"
            warning("Variables (Shot)X, (Shot)Y were renamed to point.x, point.y")
            names(x) <- dfNames
        }

        ## if AimX, AimY, AimZ -> rename to Aim.X, Aim.Y, Aim.Z
        if(all(c("aimx", "aimy") %in% dfNames)) {
            dfNames[dfNames %in% "aimx"] <- "aim.x"
            dfNames[dfNames %in% "aimy"] <- "aim.y"
            dfNames[dfNames %in% "aimz"] <- "aim.z"
            warning("Variables AimX, AimY were renamed to aim.x, aim.y")
            names(x) <- dfNames
        }

        x
    }

    DFs <- lapply(DFs, replaceXY)

    if(combine) {
        return(combineData(DFs))
    } else {
        return(DFs)
    }
}

## read files from OnTarget-output (version 1.*)
readDataOT1 <-
function(fPath=".", fNames, fPat, combine=TRUE) {
    files <- getFileNames(fPath, fNames, fPat)
    if(length(files) < 1L) { stop("No files were selected") }

    ## assumed variables: Project Title, Group, Ammunition, Distance,
    ## Aim X, Aim Y, Center X, Center Y, Point X, Point Y
    ## 10 fields + trailing tab = 11
    nFields <- unlist(lapply(files, function(x) count.fields(x, sep="\t")))
    if(!all(nFields == 11)) {
        stop(c("It appears at least one file does not contain exactly\n",
               "the required set of 10 variables - see help(readDataOT1)\n",
               "maybe you should use readDataMisc() instead"))
    }

    ## read in files into a list of data frames
    DFs <- lapply(files, function(f) {
        DF <- read.delim(f, colClasses=c("character", "factor", "character",
            "numeric", "numeric", "numeric", "numeric", "numeric",
            "numeric", "numeric", "NULL"), strip.white=TRUE)
        DF$file <- basename(tools::file_path_sans_ext(f))  # add filename
        setNames(DF, tolower(names(DF)))  # convert var names to lower case
    })

    names(DFs) <- paste("file", seq_along(DFs), sep="")  # name them

    ##  build shared set of variable names
    varNames <- Reduce(intersect, lapply(DFs, names))

    ## make sure that the data frames all have the correct variables
    wants <- c("group", "distance", "aim.x", "aim.y", "point.x", "point.y")
    has   <- wants %in% varNames
    if(!all(has)) {
        warning(c("At least one file is missing variable(s)\n",
                  paste(wants[!has], collapse= " "),
                  "\nthat may be required later by analysis functions"))
    }

    if(combine) {
        return(combineData(DFs))
    } else {
        return(DFs)
    }
}

## read files from OnTarget-output version 2.*, 3.7*, 3.8*
readDataOT2 <-
function(fPath=".", fNames, fPat, combine=TRUE) {
    files <- getFileNames(fPath, fNames, fPat)
    if(length(files) < 1L) { stop("No files were selected") }

    ## assumed variables: Project Title, Group, Ammunition, Distance,
    ## Aim X, Aim Y, Center X, Center Y, Point X, Point Y, Velocity (optional)
    ## 10 or 11 fields
    nFields <- unlist(lapply(files, function(x) count.fields(x, sep=",")))
    colClasses <- if(all(nFields == 10)) {
        c("character", "factor", "character", "numeric", "numeric",
          "numeric", "numeric", "numeric", "numeric", "numeric")
    } else if(all(nFields == 11)) {
        c("character", "factor", "character", "numeric", "numeric",
          "numeric", "numeric", "numeric", "numeric", "numeric", "numeric")
    } else {
        stop(c("It appears at least one file does not contain exactly\n",
               "the required set of variables - see help(readDataOT2)\n",
               "maybe you should use readDataMisc() instead"))
    }

    ## read in files into a list of data frames
    DFs <- lapply(files, function(f) {
        DF      <- read.csv(f, colClasses=colClasses, strip.white=TRUE)
        DF$file <- basename(tools::file_path_sans_ext(f))  # add filename
        setNames(DF, tolower(names(DF))) # convert var names to lower case
    })

    names(DFs) <- paste("file", seq_along(DFs), sep="")  # name them

    ##  build shared set of variable names
    varNames <- Reduce(intersect, lapply(DFs, names))

    ## make sure that the data frames all have the correct variables
    wants <- c("group", "distance", "aim.x", "aim.y", "point.x", "point.y")
    has   <- wants %in% varNames
    if(!all(has)) {
        warning(c("At least one file is missing variable(s)\n",
                  paste(wants[!has], collapse= " "),
                  "\nthat may be required later by analysis functions"))
    }

    if(combine) {
        return(combineData(DFs))
    } else {
        return(DFs)
    }
}
