.GET_default <-
function(url)
{
    response <- curl_fetch_memory(url)
    list(header = .status_and_headers(response),
         body = rawToChar(response$content))
}

.GET_using_gzip_compression <-
function(url)
{    
    h <- new_handle("Accept_Encoding" = "gzip")
    response <- curl_fetch_memory(url, h)
    header <- .status_and_headers(response)
    if(!grepl("gzip", header["Content-Encoding"]) ||
        inherits(tryCatch(body <-
                              memDecompress(response$content, "gzip",
                                            TRUE),
                          error = identity),
                 "error"))
        body <- rawToChar(response$content)
    list(header = header, body = body)
}

GET <-
local({
    val <- .GET_default
    function(new, gzip = FALSE) {
        if(!missing(new))
            val <<- new
        else if(!missing(gzip)) {
            if(identical(gzip, TRUE))
                val <<- .GET_using_gzip_compression
            else
                val <<- .GET_default
        } else {
            val
        }
    }
})

.status_and_headers <-
function(response)
{
    lines <- parse_headers(response$headers)
    parts <- strsplit(lines[1L], " ")[[1L]]
    s <- c("Status-Code" = parts[2L],
           "Reason-Phrase" = paste(parts[-c(1L, 2L)], collapse = " "))
    lines <- lines[-1L]
    h <- sub("[^:]+: (.*)", "\\1", lines)
    names(h) <- sub("([^:]+):.*", "\\1", lines)
    c(s, h)
}
