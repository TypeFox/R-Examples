
# A rawFormat is a list of rawBlock's

flattenBlock <- function(block) {
    UseMethod("flattenBlock")
}

flattenBlock.default <- function(block) {
    do.call("c", lapply(block, flattenBlock))
}

flattenBlock.rawBlock <- function(block) {
    list(block)
}

flattenFormat <- function(format) {
    flatFormat <- do.call("c", lapply(format, flattenBlock))
}

as.character.rawFormat <- function(x, sep1="  :  ", sep2="  |  ",
                                   blockHead=TRUE, blockChar="=", ...) {
    x$blocks <- flattenFormat(x$blocks)
    class(x) <- c("flatRawFormat", "rawFormat")
    as.character(x, sep1=sep1, sep2=sep2,
                 blockHead=blockHead, blockChar=blockChar, ...)
}

as.character.flatRawFormat <- function(x, sep1="  :  ", sep2="  |  ",
                                       blockHead=TRUE, blockChar="=", ...) {
    # Strip metadata from rawFormat to leave only a list
    # of rawBlock's
    formatOffset <- x$offset
    formatNbytes <- x$nbytes

    nblock <- length(x$blocks)

    # Generate text for each block
    blockTxt <- lapply(x$blocks, rawBlockTxt)

    # Format each block so that the separators line up
    widthChars <- getOption("width")
    offsetRange <- c(formatOffset, formatOffset + formatNbytes - 1)
    ncharOffset <- nchar(format(offsetRange)[1])
    ncharSeps <- nchar(sep1) + nchar(sep2)
    charsLeft <- widthChars - ncharOffset - ncharSeps
    # Generate initial widths for each block
    width <- vector("list", nblock)
    for (i in 1:nblock) {
        block <- x$blocks[[i]]
        blockText <- blockTxt[[i]]
        if (is.null(block$width)) {
            blockWidth <- calculateWidth(block, blockText,
                                         offsetRange, sep1, sep2)
        } else {
            blockWidth <- block$width
        }
        ncharPerMachine <- nchar(blockText$rawFiller) + 1 
        ncharPerHuman <- nchar(blockText$humanFiller) +
            nchar(blockText$humanSep)
        # For each block, store the width (num bytes)
        # and the number of chars per machine byte representation
        # and the number of chars per human byte representation
        width[[i]] <- c(blockWidth, ncharPerMachine, ncharPerHuman)
    }
    # Determine widest machine format
    maxMachineChars <- max(sapply(width, function(x) x[1]*x[2]))
    # Calculate the the space allowed for
    # human format on the line with the widest machine format
    humanChars <- charsLeft - maxMachineChars
    # Look at each block and either pad the machine format
    # or reduce the number of bytes (and pad) in order to
    # line up the sep2
    pad <- rep("", nblock)
    for (i in 1:nblock) {
        block <- x$blocks[[i]]
        blockWidth <- width[[i]]
        # Calculate the difference between the space allowed for
        # human format and the space required by this block for
        # human format
        humanDiff <- humanChars -
            (blockWidth[1] %/% block$size)*blockWidth[3]
        # If there is not enough room for the current human format
        # then we need to reduce the number of bytes
        if (humanDiff < 0) {
            # Amount of reduction
            reduction <- abs(humanDiff) %/% blockWidth[3]
            if (abs(humanDiff) %% blockWidth[3] != 0)
                reduction <- reduction + 1
            blockWidth[1] <- blockWidth[1] - reduction*block$size
            # If we reduce number of bytes to zero then
            # put it back up to 1 and issue a warning
            if (blockWidth[1] < block$size) {
                warning("Unable to align block")
                blockWidth[1] <- block$size
            }
            humanDiff <- humanChars - 
                (blockWidth[1] %/% block$size)*blockWidth[3]
        }
        # Now pad the machine format if necessary
        # humanDiff should be non-negative now!
        if (humanDiff < 0)
            stop("Houston, we have a problem ...")
        width[[i]] <- blockWidth
        machineWidth <- blockWidth[1]*blockWidth[2]
        machineDiff <- maxMachineChars - machineWidth
        if (machineDiff > 0)
            pad[i] <- paste(rep(" ", machineDiff), collapse="")
    }

    # Generate strings for each block
    formatStrings <- vector("list", nblock)
    for (i in 1:nblock) {
        block <- x$blocks[[i]]
        # Add max offset to offset sequence for formatting
        # then drop it once formatted
        offsetCol <- format(c(seq(block$offset,
                                  block$offset + block$nbytes - 1,
                                  by=width[[i]][1]),
                              offsetRange[2]))
        formatStrings[[i]] <- blockStrings(block,
                                           offsetCol[-length(offsetCol)],
                                           blockTxt[[i]],
                                           width[[i]][1],
                                           sep1, sep2, pad[i])
    }
    if (is.null(names(x$blocks))) {
        names(x$blocks) <- seq_along(x$blocks)
    }
    if (blockHead) {
        if (!is.character(blockHead))
            blockHead <- paste(rep(blockChar,
                                   ncharOffset + nchar(sep1)),
                               collapse="")
        names(formatStrings) <- paste(blockHead,
                                      names(x$blocks), sep="")
    }
    formatStrings
}

print.rawFormat <- function(x, sep1="  :  ", sep2="  |  ",
                            blockHead=TRUE, blockChar="=",
                            page=FALSE, ...) {
    rawFormat <- as.character(x, sep1=sep1, sep2=sep2,
                              blockHead=blockHead, blockChar=blockChar)
    rawFormatNames <- names(rawFormat)
    
    # View everything
    if (page) {
        tmpfile <- tempfile()
        for (i in 1:length(rawFormat)) {
            writeLines(rawFormatNames[i])
            writeLines(rawFormat[[i]], tmpfile)
        }
        file.show(tmpfile)
    } else {
        for (i in 1:length(rawFormat)) {
            cat(rawFormatNames[i], "\n")
            cat(paste(rawFormat[[i]], collapse="\n"), "\n")
        }
    }
}

