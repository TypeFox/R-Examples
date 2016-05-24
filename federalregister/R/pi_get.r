pi_get <- function(docnumber, version='v1', ...){

    # combine multiple document numbers
    docnumbers <- paste(docnumber, collapse=',')
    baseurl <- paste('https://www.federalregister.gov/api/',version,
                     '/public-inspection-documents/',docnumbers,'.json', sep='')
    
    args <- NULL
    
    r <- GET(url = paste(baseurl, args, sep=''), ...)
    stop_for_status(r)
    response <- content(r, "text")
    out <- fromJSON(response)
    if(length(strsplit(docnumbers,',')[[1]])==1){
        class(out) <- 'fedreg_document'
        out <- list(out)
    } else {
        out <- out[[2]]
        out <- lapply(out, function(x) {
            class(x) <- 'fedreg_document'
            return(x)
        })
    }
    return(out)
}
