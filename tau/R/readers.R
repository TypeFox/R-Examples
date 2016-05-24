readBytes <-
function(con)
{
    s <- rawToChar(read_all_bytes(con))
    Encoding(s) <- "bytes"
    s
}

readChars <-
function(con, encoding = "")    
{
    s <- rawToChar(read_all_bytes(con))    
    iconv(s, from = encoding, to = "UTF-8", sub = "byte")
}

read_all_bytes <-
function(con, chunksize = 2 ^ 16)
{
    if(is.character(con)) {
        return(readBin(con, raw(), file.info(con)$size))
    }

    if(!isOpen(con)) {
        open(con, "rb")
        on.exit(close(con))
    }

    bytes <- list()
    repeat {
        chunk <- readBin(con, raw(), chunksize)
        bytes <- c(bytes, list(chunk))
        if(length(chunk) < chunksize) break
    }

    unlist(bytes)
}
