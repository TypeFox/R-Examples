
# Methods for rawBlock classes (eps. print method)

# viewRaw.R does only raw blocks

# viewFormat.R does "formats" (mixture of blocks of different types)

fileSize <- function(file) {
    file.info(file)$size
}

# Store "width" and "machine" even though they are only used for printing
# because they will serve as a default "native" format for printing
# the raw block [can be overridden is as.character() and/or print()]
readRawBlock <- function(infile, width, machine,
                         type, offset, nbytes, size, endian, signed) {
    
    fileRaw <- readBin(infile, what="raw", n=nbytes)

    # get fileNum here because may require re-reading
    # the file, and we do not want to read the file in the print
    # method
    if (type == "char") {
        fileNum <- NULL
    } else {
        if (type == "int") {
            what <- "int"
        } else { # "real"
            what <- "double"
        }

        # Reset file position
        seek(infile, offset)
               
        fileNum <- readBin(infile, what=what, size=size,
                           n=nbytes %/% size,
                           endian=endian, signed=signed)
    }

    rawBlock <- list(width=width, offset=offset, nbytes=nbytes,
                     fileRaw=fileRaw, fileNum=fileNum,
                     machine=machine, type=type,
                     size=size, endian=endian, signed=signed)
    class(rawBlock) <- "rawBlock"
    rawBlock
}

# Special case for ASCII string terminated by '\n'
# (do something similar for '\0' terminated C string?)
readASCIIline <- function(infile, offset) {
    nbytes <- 0
    repeat {
        raw <- readBin(infile, what="raw", n=1)
        if (raw < 128 &&
            rawToChar(raw) == '\n')
            break
        else
            nbytes <- nbytes + 1
    }

    # Reset file location
    seek(infile, offset)
    
    readRawBlock(infile, NULL, "hex", "char", offset, nbytes + 1,
                 1, "little", TRUE)
}

blockASCII <- function(raw, nbytes, vector=TRUE) {
    # Only printing ASCII characters
    ASCIIraw <- (raw < 128 & raw > 31)[1:nbytes]
    
    ASCIItxt <- rep(".", nbytes)
    rawASCII <- rawToChar(raw[1:nbytes][ASCIIraw], multiple=TRUE)
    ASCIItxt[ASCIIraw] <- rawASCII
    if (vector)
        ASCIItxt
    else
        paste(ASCIItxt, collapse="")
}

# Extract useful "value" from block
blockValue <- function(block) {
    if (block$type == "char") {
        blockASCII(block$fileRaw, block$nbytes)
    } else { # "int" or "real"
        block$fileNum
    }
}

# Extract a null-terminated string from block
blockString <- function(block) {
    if (block$type != "char")
        stop("Invalid block type")
    # Find first null char in string
    nullPos <- match(TRUE, block$fileRaw == 0)
    # If there are no null chars return whole string
    if (is.na(nullPos)) {
        blockASCII(block$fileRaw, block$nbytes, vector=FALSE)
    } else {
        blockASCII(block$fileRaw, nullPos - 1, vector=FALSE)
    }
}

blockRaw <- function(block) {
    block$fileRaw
}

rawBlockTxt <- function(block, showHuman=TRUE) {
    with(block,
         {
             # Format the machine output
             if (machine == "binary") { 
                 fileBits <- rawToBits(fileRaw)
                 # rev() cos bits are least-significant bit first
                 # substr() cos each bit gets printed as two binary chars
                 fileRawTxt <- apply(matrix(substr(as.character(fileBits),
                                                   2 ,2),
                                            byrow=TRUE, ncol=8),
                                     1,
                                     function(x) {
                                         paste(rev(x), collapse="")
                                     })
                 rawFiller <- "        "
             } else { # "hex"
                 fileRawTxt <- as.character(fileRaw)
                 rawFiller <- "  "
             }

             # Format the human output
             if (type == "char") {
                 fileCharTxt <- blockASCII(fileRaw, nbytes)
                 
                 humanFiller <- " "
                 humanSep <- ""
                 bytesPerHuman <- 1
             } else { # int or real
                 fileCharTxt <- format(fileNum)
                 
                 humanFiller <- paste(rep(" ", nchar(fileCharTxt[1])),
                                      collapse="")
                 humanSep <- " "
                 bytesPerHuman <- size
             }
             
             # If not going to display the human column
             # "zero" the strings
             if (!showHuman) {
                 fileCharTxt <- ""
                 humanFiller <- ""
                 humanSep <- ""
             }
             
             list(fileRawTxt=fileRawTxt, rawFiller=rawFiller,
                  fileCharTxt=fileCharTxt,
                  humanFiller=humanFiller, humanSep=humanSep,
                  bytesPerHuman=bytesPerHuman)
         })
}

calculateWidth <- function(block, blockTxt, offsetRange, sep1, sep2,
                           showOffset=TRUE, showHuman=TRUE) {
    widthChars <- getOption("width")
    ncharSeps <- 0
    if (showOffset) {
        ncharOffset <- nchar(format(offsetRange)[1])
        ncharSeps <- ncharSeps + nchar(sep1)
    } else {
        ncharOffset <- 0
    }
    if (showHuman) {
        ncharSeps <- ncharSeps + nchar(sep2)
    }
    # includes " " sep
    ncharPerMachine <- nchar(blockTxt$rawFiller) + 1 
    ncharPerHuman <- nchar(blockTxt$humanFiller) + nchar(blockTxt$humanSep) 
    width <- (widthChars - ncharOffset - ncharSeps) %/%
        (ncharPerMachine*block$size + ncharPerHuman)
    if (width == 0)
        stop("Unable to fit any bytes in available width")
    width*block$size
}

blockStrings <- function(block, offsetCol, blockTxt, width,
                         sep1, sep2, pad="",
                         showOffset=TRUE, showHuman=TRUE) {
    nrow <- block$nbytes %/% width
    if (block$nbytes %% width != 0)
        nrow <- nrow + 1

    if (!showOffset) {
        sep1 <- ""
    }
    if (!showHuman) {
        sep2 <- ""
    }
    
    rawVector <- rep(blockTxt$rawFiller, nrow*width)
    rawVector[1:block$nbytes] <- blockTxt$fileRawTxt
    rawMatrix <- matrix(rawVector, byrow=TRUE, ncol=width)
             
    humanVector <- rep(blockTxt$humanFiller,
                       nrow*width %/% blockTxt$bytesPerHuman)
    humanVector[1:length(blockTxt$fileCharTxt)] <-
        blockTxt$fileCharTxt
    humanMatrix <- matrix(humanVector, byrow=TRUE,
                          ncol=width %/% blockTxt$bytesPerHuman)
    humanCol <- apply(humanMatrix, 1, paste,
                      collapse=blockTxt$humanSep)
    
    viewMatrix <- cbind(offsetCol, sep1,
                        apply(rawMatrix, 1, paste, collapse=" "),
                        pad, sep2, humanCol)
    apply(viewMatrix, 1, paste, collapse="")
}

as.character.rawBlock <- function(x, width=NULL, machine=NULL,
                                  sep1="  :  ", sep2="  |  ",
                                  showOffset=TRUE, showHuman=TRUE, ...) {
    if (!is.null(width))
        x$width <- width
    
    if (!is.null(machine))
        x$machine <- machine

    with(x,
         {
             if(!(machine %in% c("hex", "binary")))
                 stop("Invalid machine format")
    
             if ((type != "char") &&
                 !is.null(width) &&
                 (width %% size != 0))
                 stop("Incompatible width and human format")

             blockTxt <- rawBlockTxt(x, showHuman)
             
             # If width is NULL, use options("width") to determine
             # width (to fit)
             if (is.null(width)) {
                 width <- calculateWidth(x, blockTxt,
                                         c(offset, offset + nbytes - 1),
                                         sep1, sep2,
                                         showOffset, showHuman)
             }
             
             # Format the offset
             if (showOffset) {
                 offsetCol <- format(seq(offset, offset + nbytes - 1,
                                         by=width))
             } else {
                 offsetCol <- ""
             }
             
             # Put offset, machine text, and human text together
             blockStrings(x, offsetCol, blockTxt, width, sep1, sep2, 
                          showOffset=showOffset, showHuman=showHuman)
         })
}

print.rawBlock <- function(x, width=NULL, machine=NULL,
                           sep1="  :  ", sep2="  |  ",
                           showOffset=TRUE, showHuman=TRUE,
                           page=FALSE, ...) {
    rawBlock <- as.character(x, width=width, machine=machine,
                             sep1=sep1, sep2=sep2,
                             showOffset=showOffset, showHuman=showHuman)
    
    # View everything
    if (page) {
        tmpfile <- tempfile()
        writeLines(rawBlock, tmpfile)
        file.show(tmpfile)
    } else {
        cat(paste(rawBlock, collapse="\n"), "\n")
    }
}


