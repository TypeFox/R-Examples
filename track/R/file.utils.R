##  getAbsolutePath() is from Henrik Bengtsson's <henrikb@braju.com> R.utils package, Version 0.8.0
##  Copied here (with dependencies), rather than using the R.utils package,
##  because the R.utils package requires the R.oo package,
##  which redefines a whole lot of stuff, which I'd rather not do.

##  TAP: * changed [A-Z] to [A-Za-z] (want to change to [[:alpha:]]?)
##       (for detecting drive letters in paths)
##       * call dir.exists() instead of file.exists() on directories
##

extract.file.utils.funcs <- function(file=NULL) {
    # have to have R.utils attached ahead of track when running this func
    if (any(!is.element(c("package:R.utils"), search())))
        stop("must have package R.utils loaded")
    if (!is.null(file)) {
        sink(file)
        on.exit(sink())
    }
    catn <- function(x) cat(paste(x, "\n", sep=""), sep="")
    cat("\nfilePath <-\n    ")
    catn(deparse(get("filePath.default", pos="package:R.utils"), control=character(0), width.cutoff=75))
    cat("\ngetAbsolutePath <-\n    ")
    catn(deparse(get("getAbsolutePath", pos="package:R.utils"), control=character(0), width.cutoff=75))
    cat("\nisAbsolutePath <-\n    ")
    catn(deparse(get("isAbsolutePath", pos="package:R.utils"), control=character(0), width.cutoff=75))
    cat("\nreadWindowsShortcut <-\n    ")
    catn(deparse(get("readWindowsShortcut", pos="package:R.utils"), control=character(0), width.cutoff=75))
    cat("\nintToBin <-\n    ")
    catn(deparse(get("intToBin", pos="package:R.utils"), control=character(0), width.cutoff=75))
    cat("\nintToChar <- function(x) rawToChar(as.raw(x), multiple = TRUE)")
}

filePath <-
    function (..., fsep = .Platform$file.sep, removeUps = TRUE, expandLinks = c("none",
    "any", "local", "relative", "network"), mustExist = FALSE, verbose = FALSE)
{
    removeEmptyDirs <- function(pathname) {
        isOnNetworkB <- (regexpr("^\\\\\\\\", pathname) != -1)
        isOnNetworkF <- (regexpr("^//", pathname) != -1)
        pathname <- gsub("///*", "/", pathname)
        pathname <- gsub("\\\\\\\\\\\\*", "\\\\", pathname)
        if (isOnNetworkB) {
            pathname <- paste("\\\\", pathname, sep = "")
            pathname <- gsub("^\\\\\\\\\\\\*", "\\\\", pathname)
        }
        if (isOnNetworkF) {
            pathname <- paste("//", pathname, sep = "")
            pathname <- gsub("^///*", "//", pathname)
        }
        pathname
    }
    removeUpsFromPathname <- function(pathname, split = FALSE) {
        if (regexpr("^[A-Za-z]:[/\\]$", pathname) != -1)
            return(gsub("\\\\", "/", pathname))
        components <- strsplit(pathname, split = "[/\\]")[[1]]
        if (length(components) > 1) {
            components <- components[components != "."]
        }
        pos <- 2
        while (pos <= length(components)) {
            if (components[pos] == ".." && components[pos - 1] != "..") {
                if (verbose) {
                  cat("Removing: ", paste(components[c(pos - 1, pos)], collapse = ", "),
                    "\n", sep = "")
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
    isEmpty <- unlist(lapply(args, FUN = function(x) (length(x) == 0)))
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
    if (regexpr("^[A-Z]:[/\\]$", pathname) != -1)
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
        isWindowsShortcut <- (regexpr("[.](lnk|LNK)$", pathname) != -1)
        if (isWindowsShortcut) {
            lnkFile <- pathname
        }
        else {
            if (dir.exists(pathname)) {
                expandedPathname <- pathname
                isFirst <- FALSE
                next
            }
            if (isFirst) {
                isFirst <- FALSE
                if (dir.exists(paste(pathname, "", sep = fsep))) {
                  expandedPathname <- pathname
                  next
                }
            }
            lnkFile <- paste(pathname, c("lnk", "LNK"), sep = ".")
            lnkFile <- lnkFile[file.exists(lnkFile)]
            if (length(lnkFile) == 0) {
                if (verbose) {
                  msg <- paste("Failed to expand pathname '", pathname0, "'. No target found for: ",
                    pathname, sep = "")
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
            pathname <- paste(expandedPathname, lnk$relativePath, sep = fsep)
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

getAbsolutePath <-
    function (pathname, workDirectory = getwd(), expandTilde = FALSE, ...)
{
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
            name <- substring(name, attr(reg, "match.length") + 1)
        if (removeSuffix)
            name <- gsub("[.][^.]*$", "", name)
        name
    }
    if (is.null(pathname))
        pathname <- "."
    if (!isAbsolutePath(pathname)) {
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
        pathname <- filePath(pathname, removeUps = TRUE)
    }
    if (expandTilde)
        pathname <- file.path(dirname(pathname), basename(pathname))
    pathname
}

isAbsolutePath <-
    function (pathname, ...)
{
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

readWindowsShortcut <-
    function (con, verbose = FALSE, ...)
{
    readByte <- function(con, n = 1) {
        readBin(con = con, what = "integer", size = 1, n = n, signed = FALSE,
            endian = "little")
    }
    readWord <- function(con, n = 1) {
        readBin(con = con, what = "integer", size = 2, n = n, signed = FALSE,
            endian = "little")
    }
    readDWord <- function(con, n = 1) {
        readBin(con = con, what = "integer", size = 4, n = n, signed = FALSE,
            endian = "little")
    }
    readQWord <- function(con, n = 1) {
        readBin(con = con, what = "integer", size = 8, n = n, signed = FALSE,
            endian = "little")
    }
    readString <- function(con, nchars = -1, unicoded = FALSE) {
        if (nchars == -1) {
            bfr <- c()
            while ((byte <- readByte(con)) != 0) {
                bfr <- c(bfr, byte)
            }
        }
        else {
            if (unicoded)
                nchars <- 2 * nchars
            bfr <- readByte(con, n = nchars)
        }
        if (unicoded)
            bfr <- bfr[bfr != 0]
        paste(intToChar(bfr), collapse = "")
    }
    if (is.character(con)) {
        con <- file(con, open = "")
    }
    if (inherits(con, "connection")) {
        if (!isOpen(con)) {
            open(con, open = "rb")
            on.exit({
                if (inherits(con, "connection") && isOpen(con)) close(con)
            })
        }
    }
    header <- list(magic = readDWord(con), guid = readByte(con, n = 16), flags = readDWord(con),
        fileAttributes = readDWord(con), creationTime = readQWord(con), modificationTime = readQWord(con),
        lastAccessTime = readQWord(con), fileLength = readDWord(con), iconNumber = readDWord(con),
        showWndValue = readDWord(con), hotKey = readDWord(con), unknown = readDWord(con,
            n = 2))
    if (verbose) {
        cat("File header read:")
        print(header)
    }
    if (header$magic != 76) {
        stop(paste("File format error: Magic dword in header is not 0000004C (76): ",
            header$magic, sep = ""))
    }
    knownGuid <- c(1, 20, 2, 0, 0, 0, 0, 0, 192, 0, 0, 0, 0, 0, 0, 70)
    if (!all.equal(header$guid, knownGuid)) {
        stop(paste("File format error: Unknown GUID: ", paste(header$guid, collapse = ","),
            sep = ""))
    }
    flags <- intToBin(header$flags)
    flags <- rev(strsplit(flags, split = "")[[1]])
    flags <- as.logical(as.integer(flags))
    if (length(flags) > 8)
        stop(paste("File format error: Too many bits in flags in header: ",
            length(flags), sep = ""))
    flags <- c(flags, rep(FALSE, length.out = 8 - length(flags)))
    names(flags) <- c("hasShellItemIdList", "pointsToFileOrDirectory", "hasDescription",
        "hasRelativePath", "hasWorkingDirectory", "hasCommandLineArguments",
        "hasCustomIcon", "unicodedStrings")
    header$flags <- flags
    rm(flags)
    if (header$flags["pointsToFileOrDirectory"]) {
        fileAttributes <- intToBin(header$fileAttributes)
        fileAttributes <- rev(strsplit(fileAttributes, split = "")[[1]])
        fileAttributes <- as.logical(as.integer(fileAttributes))
        if (length(fileAttributes) > 13)
            stop(paste("File format error: Too many bits in flags in header: ",
                length(fileAttributes), sep = ""))
        fileAttributes <- c(fileAttributes, rep(FALSE, length.out = 13 - length(fileAttributes)))
        names(fileAttributes) <- c("isReadOnly", "isHidden", "isSystemFile",
            "isVolumeLabel", "isDirectory", "isModifiedSinceLastBackup", "isEncrypted",
            "isNormal", "isTemporary", "isSparseFile", "hasReparsePointData",
            "isCompressed", "isOffline")
        header$fileAttributes <- fileAttributes
    }
    else {
        if (!all(header$fileAttributes == 0)) {
            stop("File format error: When shortcut is not pointing to a file or a directory all file attributes should be zero.")
        }
        header$fileAttributes <- NA
    }
    if (header$fileLength < 0) {
        stop(paste("File format error: File length is negative: ", header$fileLength))
    }
    if (header$flags["hasCustomIcon"]) {
    }
    else {
        if (header$iconNumber != 0)
            stop(paste("File format error: Expected zero icon number: ", header$iconNumber))
    }
    swNames <- c("SW_HIDE", "SW_NORMAL", "SW_SHOWMINIMIZED", "SW_SHOWMAXIMIZED",
        "SW_SHOWNOACTIVATE", "SW_SHOW", "SW_MINIMIZE", "SW_SHOWMINNOACTIVE",
        "SW_SHOWNA", "SW_RESTORE", "SW_SHOWDEFAULT")
    if (header$showWndValue %in% 0:10) {
        names(header$showWndValue) <- swNames[header$showWndValue + 1]
    }
    else {
        stop(paste("File format error: showWndValue in header is out of range [0:10]: ",
            header$showWndValue))
    }
    if (!all(header$unknown == 0)) {
        stop(paste("File format error: Last 2 dwords in header are not zero: ",
            header$unknown, sep = ""))
    }
    lnk <- list(header = header)
    if (header$flags["hasShellItemIdList"]) {
        bytesToRead <- readWord(con)
        if (verbose)
            cat("bytesToRead=", bytesToRead, "\n", sep = "")
        readByte(con, n = bytesToRead)
        bytesToRead <- 0
        while (bytesToRead > 0) {
            itemLength <- readWord(con)
            if (verbose)
                cat("itemLength=", itemLength, "\n", sep = "")
            bytesToRead <- bytesToRead - 2
            item <- readByte(con, n = itemLength - 2)
            print(paste(intToChar(item), collapse = ""))
            str(item)
            bytesToRead <- bytesToRead - itemLength
        }
    }
    if (header$flags["pointsToFileOrDirectory"]) {
        fileLocationInfo <- list(length = readDWord(con), firstOffset = readDWord(con),
            flags = readDWord(con), offsetLocalVolumeInfo = readDWord(con),
            offsetBasePathname = readDWord(con), offsetNetworkVolumeInfo = readDWord(con),
            offsetRemainingPathname = readDWord(con), .offset = 7 * 4)
        if (fileLocationInfo$flags %in% 0:3) {
        }
        else {
            stop(paste("File format error: Unknown volume flag: ", fileLocationInfo$flags,
                sep = ""))
        }
        flags <- intToBin(fileLocationInfo$flags)
        flags <- rev(strsplit(flags, split = "")[[1]])
        flags <- as.logical(as.integer(flags))
        flags <- c(flags, rep(FALSE, length.out = 2 - length(flags)))
        names(flags) <- c("availableOnLocalVolume", "availableOnNetworkShare")
        fileLocationInfo$flags <- flags
        if (fileLocationInfo$flags["availableOnLocalVolume"] != TRUE) {
            "Random garbage when bit 0 is clear in volume flags"[1]
        }
        if (fileLocationInfo$flags["availableOnNetworkShare"] != TRUE) {
            "Random garbage when bit 1 is clear in volume flags"[1]
        }
        if (fileLocationInfo$firstOffset != fileLocationInfo$.offset) {
            warning("File format warning: First offset in File Location Info is not 0x1C (28): ",
                fileLocationInfo$firstOffset)
            skip <- fileLocationInfo$firstOffset - fileLocationInfo$.offset
            readBin(con, what = "integer", size = 1, n = skip)
            fileLocationInfo$.offset <- fileLocationInfo$.offset + skip
        }
        if (verbose) {
            cat("File location info:\n")
            print(fileLocationInfo)
        }
        if (fileLocationInfo$flags["availableOnLocalVolume"]) {
            skip <- fileLocationInfo$offsetLocalVolumeInfo - fileLocationInfo$.offset
            readBin(con, what = "integer", size = 1, n = skip)
            fileLocationInfo$.offset <- fileLocationInfo$.offset + skip
            table <- list(length = readDWord(con), typeOfVolume = readDWord(con),
                volumeSerialNumber = readDWord(con), offsetName = readDWord(con),
                volumeLabel = "", .offset = 4 * 4)
            if (table$typeOfVolume %in% 0:6) {
                names(table$typeOfVolume) <- c("Unknown", "No root directory",
                  "Removable", "Fixed", "Remote", "CD-ROM", "Ram drive")[table$typeOfVolume +
                  1]
            }
            else {
                stop("File format error: Unknown type of volume: ", table$typeOfVolume)
            }
            if (table$offsetName != table$.offset) {
                warning("File format warning: Offset to volume name in Local Volume Table is not 0x10 (16): ",
                  table$offsetName)
                skip <- table$offsetName - table$.offset
                readBin(con, what = "integer", size = 1, n = skip)
                table$.offset <- table$.offset + skip
            }
            table$volumeLabel <- readString(con)
            table$.offset <- table$.offset + nchar(table$volumeLabel) + 1
            if (table$.offset != table$length) {
                stop("File format error: Length of structure did not match the number of bytes read.")
            }
            fileLocationInfo$.offset <- fileLocationInfo$.offset + table$.offset
            table$length <- NULL
            table$offsetName <- NULL
            table$.offset <- NULL
            fileLocationInfo$localVolumeTable <- table
            rm(table)
            if (verbose) {
                cat("File location info / Local Volume Table:\n")
                print(fileLocationInfo$localVolumeTable)
            }
            skip <- fileLocationInfo$offsetBasePathname - fileLocationInfo$.offset
            readBin(con, what = "integer", size = 1, n = skip)
            fileLocationInfo$.offset <- fileLocationInfo$.offset + skip
            fileLocationInfo$basePathname <- readString(con)
            fileLocationInfo$.offset <- fileLocationInfo$.offset + nchar(fileLocationInfo$basePathname) +
                1
            if (verbose)
                cat("basePathname='", fileLocationInfo$basePathname, "'\n",
                  sep = "")
        }
        if (fileLocationInfo$flags["availableOnNetworkShare"]) {
            skip <- fileLocationInfo$offsetNetworkVolumeInfo - fileLocationInfo$.offset
            readBin(con, what = "integer", size = 1, n = skip)
            fileLocationInfo$.offset <- fileLocationInfo$.offset + skip
            table <- list(length = readDWord(con), unknown1 = readDWord(con),
                offsetName = readDWord(con), unknown2 = readDWord(con), unknown3 = readDWord(con),
                networkShareName = "", .offset = 5 * 4)
            if (table$offsetName != table$.offset) {
                warning("File format warning: Offset to network share name in Network Volume Table is not 0x14 (20): ",
                  table$offsetName)
                readBin(con, what = "integer", size = 1, n = table$offsetName -
                  table$.offset)
            }
            table$networkShareName <- readString(con)
            table$.offset <- table$.offset + nchar(table$networkShareName) +
                1
            if (verbose) {
                cat("File location info / Network Volume Table:\n")
                print(table)
            }
            if (table$.offset != table$unknown2) {
                warning("File format error: Length of table structure did not match the number of bytes read.")
            }
            fileLocationInfo$.offset <- fileLocationInfo$.offset + table$.offset
            table$length <- NULL
            table$offsetName <- NULL
            table$unknown1 <- table$unknown2 <- table$unknown3 <- NULL
            table$.offset <- NULL
            fileLocationInfo$networkVolumeTable <- table
            rm(table)
            if (verbose) {
                cat("File location info / Network Volume Table:\n")
                print(fileLocationInfo$networkVolumeTable)
            }
        }
        skip <- fileLocationInfo$offsetRemainingPathname - fileLocationInfo$.offset
        readBin(con, what = "integer", size = 1, n = skip)
        fileLocationInfo$.offset <- fileLocationInfo$.offset + skip
        fileLocationInfo$remainingPathname <- readString(con)
        fileLocationInfo$.offset <- fileLocationInfo$.offset + nchar(fileLocationInfo$remainingPathname) +
            1
        if (fileLocationInfo$length != fileLocationInfo$.offset) {
            stop("File format error: Expected to read ", fileLocationInfo$length,
                " bytes in File Location Info structure, but read ", fileLocationInfo$.offset)
        }
        fileLocationInfo$length <- NULL
        fileLocationInfo$firstOffset <- NULL
        fileLocationInfo$offsetBasePathname <- NULL
        fileLocationInfo$offsetLocalVolumeInfo <- NULL
        fileLocationInfo$offsetNetworkVolumeInfo <- NULL
        fileLocationInfo$offsetRemainingPathname <- NULL
        fileLocationInfo$.offset <- NULL
        lnk$fileLocationInfo <- fileLocationInfo
        rm(fileLocationInfo)
    }
    else {
        lnk$fileLocationInfo <- NA
    }
    unicoded <- header$flags["unicodedStrings"]
    if (header$flags["hasDescription"]) {
        nchars <- readWord(con)
        lnk$description <- readString(con, nchars = nchars, unicoded = unicoded)
    }
    if (header$flags["hasRelativePath"]) {
        nchars <- readWord(con)
        lnk$relativePath <- readString(con, nchars = nchars, unicoded = unicoded)
    }
    if (header$flags["hasWorkingDirectory"]) {
        nchars <- readWord(con)
        lnk$workingDirectory <- readString(con, nchars = nchars, unicoded = unicoded)
    }
    if (header$flags["hasCommandLineArguments"]) {
        nchars <- readWord(con)
        lnk$commandLineArguments <- readString(con, nchars = nchars, unicoded = unicoded)
    }
    if (header$flags["hasCustomIcon"]) {
        nchars <- readWord(con)
        lnk$iconFilename <- readString(con, nchars = nchars, unicoded = unicoded)
    }
    if (header$flags["pointsToFileOrDirectory"]) {
        if (lnk$fileLocationInfo$flags["availableOnLocalVolume"]) {
            lnk$pathname <- paste(lnk$fileLocationInfo$basePathname, lnk$fileLocationInfo$remainingPathname,
                sep = "")
        }
        if (lnk$fileLocationInfo$flags["availableOnNetworkShare"]) {
            lnk$networkPathname <- paste(lnk$fileLocationInfo$networkVolumeTable$networkShareName,
                "\\", lnk$fileLocationInfo$remainingPathname, sep = "")
        }
    }
    lnk
}

intToBin <-
    function (x)
{
    y <- as.integer(x)
    class(y) <- "binmode"
    y <- as.character(y)
    dim(y) <- dim(x)
    y
}

# Re prior use of "ASCII" table
# BDR wrote on 2008-04-24 in a private message to maintainers of packages containing "\0":
# ... strings like "\0" or "\000".  In all released versions of R that has
# (AFAIK silently) parsed to "", but R 2.8.0 will give a warning.
# If you really feel you need intToChar, use
intToChar <- function(x) rawToChar(as.raw(x), multiple = TRUE)
# But it's dangerous in these days of multibyte locales -- we do have intToUtf8 + iconv.
# old code: intToChar <- function (i, ...) {ASCII[i%%256 + 1]}
