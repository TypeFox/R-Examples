PPT.getAbsolutePath<-function (pathname, workDirectory = getwd(), expandTilde = FALSE,...){




PPT.filePath<-function (..., fsep = .Platform$file.sep, removeUps = TRUE, expandLinks = c("none", 
    "any", "local", "relative", "network"), mustExist = FALSE, 
    verbose = FALSE) 
{
    removeEmptyDirs <- function(pathname) {
        isOnNetworkBwd <- (regexpr("^\\\\\\\\", pathname) != 
            -1)
        isOnNetworkFwd <- (regexpr("^//", pathname) != -1)
        pathname <- gsub("///*", "/", pathname)
        pathname <- gsub("\\\\\\\\\\\\*", "\\\\", pathname)
        if (isOnNetworkBwd) {
            pathname <- paste("\\\\", pathname, sep = "")
            pathname <- gsub("^\\\\\\\\\\\\*", "\\\\\\\\", pathname)
        }
        if (isOnNetworkFwd) {
            pathname <- paste("//", pathname, sep = "")
            pathname <- gsub("^///*", "//", pathname)
        }
        pathname
    }
    removeUpsFromPathname <- function(pathname, split = FALSE) {
        if (regexpr("^[ABCDEFGHIJKLMNOPQRSTUVWXYZ]:[/\\]$", pathname) != 
            -1) 
            return(gsub("\\\\", "/", pathname))
        components <- strsplit(pathname, split = "[/\\]")[[1]]
        if (length(components) > 1) {
            components <- components[components != "."]
        }
        pos <- 2
        while (pos <= length(components)) {
            if (components[pos] == ".." && components[pos - 1] != 
                "..") {
                if (verbose) {
                  cat("Removing: ", paste(components[c(pos - 
                    1, pos)], collapse = ", "), "\n", sep = "")
                }
                components <- components[-c(pos - 1, pos)]
                pos <- pos - 1
            }
            else {
                pos <- pos + 1
            }
        }
        if (split) {
            components
        }
        else {
            paste(components, collapse = fsep)
        }
    }
    args <- list(...)
    isEmpty <- unlist(lapply(args, FUN = function(x) (length(x) == 
        0)))
    args <- args[!isEmpty]
    args <- lapply(args, FUN = as.character)
    expandLinks <- match.arg(expandLinks)
    if (length(args) == 0) 
        return(NULL)
    pathname <- paste(args, collapse = fsep)
    pathname <- removeEmptyDirs(pathname)
    if (expandLinks == "none") {
        if (removeUps) 
            pathname <- removeUpsFromPathname(pathname)
        return(pathname)
    }
    if (regexpr("^[ABCDEFGHIJKLMNOPQRSTUVWXYZ]:[/\\]$", pathname) != 
        -1) 
        return(gsub("\\\\", "/", pathname))
    pathname0 <- pathname
    components <- removeUpsFromPathname(pathname, split = TRUE)
    isFirst <- TRUE
    expandedPathname <- NULL
    ready <- FALSE
    while (!ready) {
        if (length(components) == 0) {
            ready <- TRUE
            break
        }
        component <- components[1]
        components <- components[-1]
        if (isFirst) {
            pathname <- component
        }
        else {
            pathname <- paste(expandedPathname, component, sep = fsep)
        }
        if (verbose) 
            print(pathname)
        isWindowsShortcut <- (regexpr("[.](lnk|LNK)$", pathname) != 
            -1)
        if (isWindowsShortcut) {
            lnkFile <- pathname
        }
        else {
            if (file.exists(pathname)) {
                expandedPathname <- pathname
                isFirst <- FALSE
                next
            }
            if (isFirst) {
                isFirst <- FALSE
                if (file.exists(paste(pathname, "", sep = fsep))) {
                  expandedPathname <- pathname
                  next
                }
            }
            lnkFile <- paste(pathname, c("lnk", "LNK"), sep = ".")
            lnkFile <- lnkFile[file.exists(lnkFile)]
            if (length(lnkFile) == 0) {
                if (verbose) {
                  msg <- paste("Failed to expand pathname '", 
                    pathname0, "'. No target found for: ", pathname, 
                    sep = "")
                  cat(msg, "\n")
                }
                break
            }
            lnkFile <- lnkFile[1]
        }
        tryCatch({
            lnk <- readWindowsShortcut(lnkFile)
        }, error = function(ex) {
            if (verbose) {
                msg <- paste("Invalid Windows shortcut found when expanding pathname '", 
                  pathname0, "': ", lnkFile, sep = "")
                cat(msg, "\n")
                print(ex)
            }
            ready <<- TRUE
        })
        if (ready) 
            break
        pathname <- NULL
        if (expandLinks == "any") {
            if (lnk$fileLocationInfo$flags["availableOnLocalVolume"]) {
                pathname <- lnk$pathname
                if (is.null(pathname)) 
                  pathname <- lnk$relativePath
            }
            if (lnk$fileLocationInfo$flags["availableOnNetworkShare"]) {
                if (is.null(pathname)) 
                  pathname <- lnk$networkPathname
            }
        }
        else if (expandLinks == "local") {
            if (lnk$fileLocationInfo$flags["availableOnLocalVolume"]) {
                pathname <- lnk$pathname
            }
        }
        else if (expandLinks %in% c("relative")) {
            if (is.null(expandedPathname)) 
                expandedPathname <- removeUpsFromPathname(pathname0)
            pathname <- paste(expandedPathname, lnk$relativePath, 
                sep = fsep)
            if (removeUps) 
                pathname <- removeUpsFromPathname(pathname)
        }
        else if (expandLinks %in% c("network")) {
            if (lnk$fileLocationInfo$flags["availableOnNetworkShare"]) {
                pathname <- lnk$networkPathname
            }
        }
        if (is.null(pathname)) {
            if (verbose) {
                msg <- paste("No target found in Windows shortcut when expanding pathname '", 
                  pathname0, "': ", lnkFile, sep = "")
                cat(msg, "\n")
            }
            break
        }
        expandedPathname <- pathname
    }
    if (length(components) > 0) {
        if (mustExist) {
            pathname <- pathname0
        }
        else {
            pathname <- paste(pathname, paste(components, collapse = fsep), 
                sep = fsep)
        }
    }
    if (removeUps) 
        pathname <- removeUpsFromPathname(pathname)
    pathname
}


PPT.isAbsolutePath<-function (pathname, ...){
    if (length(pathname) == 0) 
        return(FALSE)
    pathname <- as.character(pathname)
    if (regexpr("^~", pathname) != -1) 
        return(TRUE)
    if (regexpr("^.:(/|\\\\)", pathname) != -1) 
        return(TRUE)
    components <- strsplit(pathname, split = "[/\\]")[[1]]
    if (length(components) == 0) 
        return(FALSE)
    (components[1] == "")
}




    getName <- function(pathname, removeSuffix = FALSE, ...) {
        components <- strsplit(pathname, split = "[/\\]")[[1]]
        len <- length(components)
        if (len == 0) 
            return("")
        name <- components[len]
        if (name == ".") 
            return("")
        reg <- regexpr("^.:", name)
        if (reg != -1) 
            name <- substring(name, attr(reg, "match.length") + 
                1)
        if (removeSuffix) 
            name <- gsub("[.][^.]*$", "", name)
        name
    }
    if (is.null(pathname)) 
        pathname <- "."
    if (!PPT.isAbsolutePath(pathname)) {
        workDirectory <- strsplit(workDirectory, split = "[/\\]")[[1]]
        name <- getName(pathname)
        if (name == "" || name == ".") 
            name <- NULL
        pathname <- strsplit(pathname, split = "[/\\]")[[1]]
        len <- length(pathname)
        if (len != 0) 
            pathname <- pathname[-len]
        pathname <- c(workDirectory, pathname, name)
        pathname <- paste(pathname, sep = "", collapse = .Platform$file.sep)
        pathname <- PPT.filePath(pathname, removeUps = TRUE)
    }
    if (expandTilde) 
        pathname <- file.path(dirname(pathname), basename(pathname))
    #pathname <- gsub("//*", "/", pathname) #Commented out to deal with network drives.
    pathname
}
