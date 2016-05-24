readWave <- 
function(filename, from = 1, to = Inf, 
    units = c("samples", "seconds", "minutes", "hours"), header = FALSE, toWaveMC = NULL){
        
    if(!is.character(filename))
        stop("'filename' must be of type character.")
    if(length(filename) != 1)
        stop("Please specify exactly one 'filename'.")
    if(!file.exists(filename))
        stop("File '", filename, "' does not exist.")
    if(file.access(filename, 4))
        stop("No read permission for file ", filename)

    ## Open connection
    con <- file(filename, "rb")
    on.exit(close(con)) # be careful ...
    int <- integer()
    
    ## Reading in the header:
    RIFF <- readChar(con, 4)
    file.length <- readBin(con, int, n = 1, size = 4, endian = "little")
    WAVE <- readChar(con, 4)
    if (!(RIFF == "RIFF" && WAVE == "WAVE"))
        warning("Looks like '", filename, "' is not a valid wave file.")
   
    FMT <- readChar(con, 4)    
    bext <- NULL
    ## extract possible bext information, if header = TRUE
    if (header && (tolower(FMT) == "bext")){
        bext.length <- readBin(con, int, n = 1, size = 4, endian = "little")
        bext <- sapply(seq(bext.length), function(x) readChar(con, 1, useBytes=TRUE))
        bext[bext==""] <- " "
        bext <- paste(bext, collapse="")
        FMT <- readChar(con, 4)
    }
        
    ## waiting for the fmt chunk
    i <- 0
    while(FMT != "fmt "){
        i <- i+1
        belength <- readBin(con, int, n = 1, size = 4, endian = "little")
        seek(con, where = belength, origin = "current")
        FMT <- readChar(con, 4)
        if(i > 5) stop("There seems to be no 'fmt ' chunk in this Wave (?) file.")
    }
    fmt.length <- readBin(con, int, n = 1, size = 4, endian = "little")
    pcm <- readBin(con, int, n = 1, size = 2, endian = "little", signed = FALSE)
    ## FormatTag: only WAVE_FORMAT_PCM (0,1), WAVE_FORMAT_IEEE_FLOAT (3), WAVE_FORMAT_EXTENSIBLE (65534, determined by SubFormat)
    if(!(pcm %in% c(0, 1, 3, 65534)))
        stop("Only uncompressed PCM and IEEE_FLOAT Wave formats supported")
    channels <- readBin(con, int, n = 1, size = 2, endian = "little")
    sample.rate <- readBin(con, int, n = 1, size = 4, endian = "little")
    bytes.second <- readBin(con, int, n = 1, size = 4, endian = "little")
    block.align <- readBin(con, int, n = 1, size = 2, endian = "little")
    bits <- readBin(con, int, n = 1, size = 2, endian = "little")
    if(!(bits %in% c(8, 16, 24, 32, 64)))
        stop("Only 8-, 16-, 24-, 32- or 64-bit Wave formats supported")
    ## non-PCM (chunk size 18 or 40) 
    if(fmt.length >= 18)    
        cbSize <- readBin(con, int, n = 1, size = 2, endian = "little")
    ## chunk size 40 (extension 22)
    if(exists("cbSize") && cbSize == 22 && fmt.length == 40){
        validBits <- readBin(con, int, n = 1, size = 2, endian = "little")
        dwChannelMask <- readBin(con, int, n = 1, size = 4, endian = "little")    
        channelNames <- MCnames[as.logical(intToBits(dwChannelMask)),"name"]
        SubFormat <- readBin(con, int, n = 1, size = 2, endian = "little", signed = FALSE)
        x <- readBin(con, "raw", n=14)
    }
    if(exists("SubFormat") && !(SubFormat %in% c(0, 1, 3)))
        stop("Only uncompressed PCM and IEEE_FLOAT Wave formats supported")
    #if(fmt.length > 26)
    #    seek(con, where = fmt.length - 26, origin = "current")
    
    ## fact chunk
#    if((pcm %in% c(0, 3)) || (pcm = 65534 && SubFormat %in% c(0, 3))) {
#      fact <- readChar(con, 4)
#      fact.length <- readBin(con, int, n = 1, size = 4, endian = "little")
#      dwSampleLength <- readBin(con, int, n = 1, size = 4, endian = "little")
#    }
    
    DATA <- readChar(con, 4)
    ## waiting for the data chunk    
    i <- 0    
    while(DATA != "data"){
        i <- i+1
        belength <- readBin(con, int, n = 1, size = 4, endian = "little")
        seek(con, where = belength, origin = "current")
        DATA <- readChar(con, 4)
        if(i > 5) stop("There seems to be no 'data' chunk in this Wave (?) file.")
    }
    data.length <- readBin(con, int, n = 1, size = 4, endian = "little")
    bytes <- bits/8
    if(((sample.rate * block.align) != bytes.second) || 
        ((channels * bytes) != block.align))
            warning("Wave file '", filename, "' seems to be corrupted.")

    ## If only header required: construct and return it
    if(header){
        return(c(list(sample.rate = sample.rate, channels = channels, 
            bits = bits, samples = data.length / (channels * bytes)),
            if(!is.null(bext)) list(bext = bext)))
    }

    ## convert times to sample numbers
    fctr <- switch(match.arg(units),
                   samples = 1,
                   seconds = sample.rate,
                   minutes = sample.rate * 60,
                   hours = sample.rate * 3600)
    if(fctr > 1) {
        from <- round(from * fctr + 1)
        to <- round(to * fctr)
    } 

    ## calculating from/to for reading in sample data    
    N <- data.length / bytes
    N <- min(N, to*channels) - (from*channels+1-channels) + 1
    seek(con, where = (from - 1) * bytes * channels, origin = "current")

    ## reading in sample data
    ## IEEE FLOAT 
    if(pcm == 3 || (exists("SubFormat") && SubFormat==3)){
                sample.data <- readBin(con, "numeric", n = N, size = bytes, 
                    endian = "little")        
    } else {
        ## special case of 24 bits
        if(bits == 24){
            sample.data <- readBin(con, int, n = N * bytes, size = 1, 
                signed = FALSE, endian = "little")
            sample.data <- as.vector(t(matrix(sample.data, nrow = 3)) %*% 256^(0:2))
            sample.data <- sample.data - 2^24 * (sample.data >= 2^23)
        } else {
            sample.data <- readBin(con, int, n = N, size = bytes, 
                signed = (bytes != 1), endian = "little")
        }
    }
    
    ## output to WaveMC if selected by the user or if dwChannelMask suggests this is a multichannel Wave
    toWaveMC <- if(pcm != 65534 || (exists("dwChannelMask") && dwChannelMask %in% c(1,3))) isTRUE(toWaveMC) else TRUE  
    
    if(toWaveMC){
        ## Constructing the WaveMC object: 
        object <- new("WaveMC", samp.rate = sample.rate, bit = bits, 
            pcm = !(pcm == 3 || (exists("SubFormat") && SubFormat==3)))
        object@.Data <- matrix(sample.data, ncol = channels, byrow=TRUE)
        if(exists("channelNames")) {
            if((lcN <- length(channelNames)) < channels)
                channelNames <- c(channelNames, paste("C", (lcN+1):channels, sep=""))
            colnames(object@.Data) <- channelNames
        }
    } else {
        ## Constructing the Wave object: 
        object <- new("Wave", stereo = (channels == 2), samp.rate = sample.rate, bit = bits, 
            pcm = !(pcm == 3 || (exists("SubFormat") && SubFormat==3)))
        if(channels == 2) {
            sample.data <- matrix(sample.data, nrow = 2)
            object@left <- sample.data[1, ]
            object@right <- sample.data[2, ]
        } else {
            object@left <- sample.data
        }
    }
    
    ## Return the Wave object
    return(object)
}
