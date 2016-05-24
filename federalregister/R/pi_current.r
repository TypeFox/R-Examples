pi_current <- function(version='v1', ...){

    baseurl <- paste('https://www.federalregister.gov/api/',version,
                     '/public-inspection-documents/current.json', sep='')
    args <- NULL
    
    r <- GET(url = paste(baseurl, args, sep=''), ...)
    stop_for_status(r)
    response <- content(r, "text")
    out <- fromJSON(response)
    out <- lapply(out[[2]], function(x) {
        class(x) <- 'fedreg_document'
        return(x)
    })
    return(out)
}
