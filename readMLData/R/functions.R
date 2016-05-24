prepareDSList <- function(pathData, pathDescription)
{
    pathData <- path.expand(pathData)
    pathDescription <- path.expand(pathDescription)
    stopifnot(file.info(pathData)$isdir)
    stopifnot(file.info(pathDescription)$isdir)
    dsList <- readDSListFromXML(paste(pathDescription, "contents.xml", sep="/"))
    dsList <- data.frame(pathData, pathDescription, dsList)
    dsList$available <- getAvailable(dsList, seq.int(length.out=nrow(dsList)), asLogical=TRUE)
    OK <- checkConsistency(dsList)
    if (OK) {
        cat(sum(dsList$available), "data sets available\n")
    }
    dsList
}

dsSearch <- function(dsList, id, searchField=c("identification", "fullName", "dirName", "files"),
            searchType=c("exact", "prefix", "suffix", "anywhere"), caseSensitive=FALSE)
{
    stopifnot(is.data.frame(dsList))
    if (is.numeric(id)) {
        stopifnot(all(id >= 1))
        stopifnot(all(id <= nrow(dsList)))
        stopifnot(all(id == floor(id)))
        ind <- id
        if (length(ind) == 0) {
            cat("no data set found for", id, "\n")
            return(invisible(NULL))
        }
        out <- data.frame(ind=ind, identification=dsList$identification[ind], stringsAsFactors=FALSE)
    } else {
        stopifnot(is.character(id))
        stopifnot(length(id) == 1)
        searchField <- match.arg(searchField)
        searchType <- match.arg(searchType)
        field <- dsList[[searchField]]
        if (!caseSensitive) {
            field <- tolower(field)
            id <- tolower(id)
        }
        if (searchType == "exact") {
            ind <- which(field == id)
        } else if (searchType == "prefix") {
            field <- substr(field, 1, nchar(id))
            ind <- which(field == id)
        } else if (searchType == "suffix") {
            field <- substr(field, nchar(field) - nchar(id) + 1, nchar(field))
            ind <- which(field == id)
        } else {
            ind <- grep(id, field, fixed=TRUE)
        }
        if (length(ind) == 0) {
            cat("no data set found for", id, "\n")
            return(invisible(NULL))
        }
        out <- data.frame(ind=ind, identification=dsList$identification[ind], stringsAsFactors=FALSE)
        if (searchField == "fullName") {
            out <- data.frame(out, fullName=dsList$fullName[ind])
        } else if (searchField == "dirName") {
            out <- data.frame(out, dirName=dsList$dirName[ind])
        } else if (searchField == "files") {
            out <- data.frame(out, dirName=dsList$files[ind])
        }
    }
    out
}

getData <- function(x)
{
    stop("a dummy function, which should be masked by a sourced specific function")
}

dsRead <- function(dsList, id, responseName=NULL, originalNames=TRUE, deleteUnused=TRUE, keepContents=FALSE)
{
    stopifnot(is.data.frame(dsList))
    id <- getIndex(dsList, id)
    identification <- dsList$identification[id]
    functionsFile <- paste(dsList$pathDescription[id], "/functions.R", sep="")
    if (file.access(functionsFile) == 0) {
        source(functionsFile, local=TRUE)
    }
    source(paste(dsList$pathDescription[id], "/scripts/", identification, ".R", sep=""), local=TRUE)
    result <- try( dat <- getData(paste(dsList$pathData[id], dsList$dirName[id], sep="/")) )
    if (inherits(result, "try-error")) {
        cat("An error occured, when reading the data.\n")
        return(invisible(NULL))
    }
    rm(getData)
    row.names(dat) <- 1:nrow(dat)
    if (ncol(dat) != dsList$originalColsNumber[id]) {
        stop("inconsistent originalColsNumber")
    }
    if (nrow(dat) != dsList$cases[id]) {
        stop("inconsistent number of cases")
    }
    out.cols <- seq.int(length.out=ncol(dat))
    out.types <- split.comma(dsList$originalColsType[id])
    if (originalNames && dsList$originalColsNames[id] != "") {
        out.names <- split.comma(dsList$originalColsNames[id])
    } else {
        out.names <- paste("V", out.cols, sep="")
    }
    if (!is.null(responseName)) {
        responsePos <- dsList$responsePos[id]
        out.names[responsePos] <- responseName
    }
    names(dat) <- out.names
    if (deleteUnused && !keepContents)
    {
        delete <- as.integer(split.comma(dsList$delete[id]))
        if (length(delete) != 0) {
            out.cols <- out.cols[ - delete]
            out.types <- out.types[out.cols]
            out.names <- out.names[out.cols]
            dat <- dat[out.cols]
        }
    }
    if (!keepContents) {
        for (j in seq.int(length.out=ncol(dat))) {
            if (out.types[j] != "n")
                dat[,j] <- factor(dat[,j])
        }
    }
    dat
}

getAvailable <- function(dsList, id=NULL, asLogical=FALSE)
{
    stopifnot(is.data.frame(dsList))
    if (is.null(id)) {
        id <- seq.int(length.out=nrow(dsList))
    }
    if (is.character(id)) {
        id <- match(id, dsList$identification)
    }
    out <- rep(NA, times=length(id))
    for (i in id) {
        files <- split.comma(dsList$files[i])
        files <- paste(dsList$pathData[i], dsList$dirName[i], files, sep="/")
        out[i] <- all(file.access(files) == 0)
    }
    if (asLogical) {
        out
    } else {
        dsList$identification[id][out]
    }
}

runCommand <- function(dsList, id, command, fileName)
{
    setwd(paste(dsList$pathData[id], dsList$dirName[id], sep="/"))
    links <- readLines(fileName)
    isComment <- grepl("^ *#", links)
    address <- links[!isComment]
    file <- basename(address)
    doDownload <- file.access(file) == -1
    cat("\n")
    for (i in seq.int(along.with=address)) {
        if (doDownload[i]) {
            system(paste(command, address[i]))
            cat("\n")
        } else {
            cat("File", file[i], "already exists\n")
            cat("\n")
        }
    }
    if (any(isComment)) {
        cat("The data set has the following downloading comments attached.\n")
        cat("\n")
        cat(paste(links[isComment], "\n", sep=""))
        cat("\n")
    }
}

dsDownload <- function(dsList, id, command, fileName)
{
    stopifnot(is.data.frame(dsList))
    id <- getIndex(dsList, id)
    if (is.null(id)) return(invisible(NULL))
    current <- getwd()
    out <- try( runCommand(dsList, id, command, fileName) )
    setwd(current)
}

