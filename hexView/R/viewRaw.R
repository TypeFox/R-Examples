
readRaw <- function(file, width=NULL,
                    offset=0, nbytes=NULL,
                    machine="hex", human="char",
                    # R defaults
                    size=switch(human, char=1, int=4, real=8),
                    # These two only used if human != "char"
                    endian=.Platform$endian, signed=TRUE) {
    fileSize <- fileSize(file)

    infile <- file(file, "rb")
    on.exit(close(infile))

    if (offset < 0 || offset > fileSize)
        stop("Invalid offset")
    if (!is.null(nbytes) &&
        (nbytes <= 0 || offset + nbytes > fileSize))
        stop("Invalid number of bytes")
    
    if (offset > 0)
        seek(infile, offset)

    if (is.null(nbytes))
        nbytes <- fileSize - offset

    if(!(machine %in% c("hex", "binary")))
        stop("Invalid machine format")
    
    if(!(human %in% c("char", "int", "real")))
        stop("Invalid human format")

    readRawBlock(infile, width, machine, human, offset, nbytes,
                 size, endian, signed)
}

viewRaw <- function(..., page=FALSE) {
    print(readRaw(...), page=page)
}
