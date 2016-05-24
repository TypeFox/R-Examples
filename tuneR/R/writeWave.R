writeWave <- 
function(object, filename, extensible = TRUE) {
    if(!is(object, "WaveGeneral")) 
        stop("'object' needs to be of class 'Wave' or 'WaveMC'")
    validObject(object)
    if(is(object, "Wave")){
        object <- as(object, "WaveMC")
        colnames(object) <- c("FL", if(ncol(object) > 1) "FR")
    }
    if(ncol(object) > 2 && !extensible)
        stop("Objects with more than two columns (multi channel) can only be written to a Wave extensible format file.")

    cn <- colnames(object)
    if((length(cn) != ncol(object) || !all(cn %in% MCnames[["name"]])) || any(duplicated(cn)))
        stop("colnames(object) must be specified and must uniquely identify the channel ordering for WaveMC objects, see ?MCnames for possible channels")
    cnamesnum <- as.numeric(factor(colnames(object), levels=MCnames[["name"]]))
    if(is.unsorted(cnamesnum))
        object <- object[,order(cnamesnum)]
    dwChannelMask <- sum(2^(cnamesnum - 1))  ##

    l <- length(object)
    sample.data <- t(object@.Data)
    dim(sample.data) <- NULL
    
    ## PCM or IEEE FLOAT
    pcm <- object@pcm                                 
    
    if(pcm) {                                                                                     
      if((object@bit == 8) && ( (max(sample.data) > 255) || (min(sample.data) < 0) ))              
          stop("for 8-bit Wave files, data range is supposed to be in [0, 255], see ?normalize")
      if((object@bit == 16) && ( (max(sample.data) > 32767) || (min(sample.data) < -32768)))
          stop("for 16-bit Wave files, data range is supposed to be in [-32768, 32767], see ?normalize")
      if((object@bit == 24) && ( (max(sample.data) > 8388607) || (min(sample.data) < -8388608)))
          stop("for 24-bit Wave files, data range is supposed to be in [-8388608, 8388607], see ?normalize")
      if((object@bit == 32) && ( (max(sample.data) > 2147483647) || (min(sample.data) < -2147483648)))
          stop("for 32-bit Wave files, data range is supposed to be in [-2147483648, 2147483647], see ?normalize")
      if(any(sample.data %% 1 > 0))
          warning("channels' data will be rounded to integers for writing the wave file")
    } else {                                                                                    
      if( (max(sample.data) > 1) || (min(sample.data) < -1) )
          stop("for IEEE float Wave files, data range is supposed to be in [-1,1], see ?normalize") 
    }
    
    # Open connection
    con <- file(filename, "wb")
    on.exit(close(con)) # be careful ...
        
    # Some calculations:
    byte <- as.integer(object@bit / 8)
    channels <- ncol(object)
    block.align <- channels * byte
    bytes <- l * byte * channels
            
    ## Writing the header:
    # RIFF
    writeChar("RIFF", con, 4, eos = NULL) 
    writeBin(as.integer(bytes + if(extensible) 72 else 36), con, size = 4, endian = "little") # cksize RIFF
    # WAVE
    writeChar("WAVE", con, 4, eos = NULL)
    # fmt chunk
    writeChar("fmt ", con, 4, eos = NULL)
    if(extensible) { # cksize format chunk
        writeBin(as.integer(40), con, size = 4, endian = "little") 
    } else {
        writeBin(as.integer(16), con, size = 4, endian = "little")
    }    
    if(!extensible) { # wFormatTag
        writeBin(as.integer(if(pcm) 1 else 3), con, size = 2, endian = "little")
    } else {
        writeBin(as.integer(65534), con, size = 2, endian = "little") # wFormatTag: extensible   
    }
    writeBin(as.integer(channels), con, size = 2, endian = "little") # nChannels
    writeBin(as.integer(object@samp.rate), con, size = 4, endian = "little") # nSamplesPerSec
    writeBin(as.integer(object@samp.rate * block.align), con, size = 4, endian = "little") # nAvgBytesPerSec
    writeBin(as.integer(block.align), con, size = 2, endian = "little") # nBlockAlign
    writeBin(as.integer(object@bit), con, size = 2, endian = "little") # wBitsPerSample
    # extensible
    if(extensible) {
        writeBin(as.integer(22), con, size = 2, endian = "little") # cbsize extensible
        writeBin(as.integer(object@bit), con, size = 2, endian = "little") # ValidBitsPerSample
        writeBin(as.integer(dwChannelMask), con, size = 4, endian = "little") #  dbChannelMask
        writeBin(as.integer(if(pcm) 1 else 3), con, size = 2, endian = "little") # SubFormat 1-2
        writeBin(as.raw(c(0,   0,   0,  0,  16,   0, 128,   0 ,  0, 170,   0,  56, 155, 113)), con) # SubFormat 3-16
        # fact
        writeChar("fact", con, 4, eos = NULL)
        writeBin(as.integer(4), con, size = 4, endian = "little") # cksize fact chunk
        writeBin(as.integer(l), con, size = 4, endian = "little") # dwSampleLength
    }
    # data
    writeChar("data", con, 4, eos = NULL)
    writeBin(as.integer(bytes), con, size = 4, endian = "little")

    # Write data:
    # PCM format
    if(pcm) { 
      if(byte == 3){
          sample.data <- sample.data + 2^24 * (sample.data < 0)
          temp <- sample.data %% (256^2)
          sample.data <- sample.data %/% 256^2
          a2 <- temp %/% 256
          temp <- temp %%  256
          writeBin(as.integer(rbind(temp, a2, sample.data)), con, size = 1, endian = "little")
      } else {
          writeBin(as.integer(sample.data), con, size = byte, endian = "little")
      }
    } else {
      writeBin(as.numeric(sample.data), con, size = byte, endian = "little")
    }

    invisible(NULL)
}
