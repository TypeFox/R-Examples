fr_get <- function(docnumber, fields=NULL, version='v1', ...){

    # combine multiple document numbers
    docnumbers <- paste(docnumber, collapse=',')
    baseurl <- paste('https://www.federalregister.gov/api/',version,'/articles/',docnumbers,'.json', sep='')
    
    if(!is.null(fields))
        fields <- curl_escape(paste('fields[]',fields,sep='=', collapse='&'))
    
    args <- fields
    
    
    r <- GET(url = paste(baseurl, args, sep='?'), ...)
    stop_for_status(r)
    response <- content(r, "text")
    out <- fromJSON(response)
    if(length(strsplit(docnumbers,',')[[1]])==1){
        class(out) <- 'fedreg_document'
        out <- list(out)
    } else {
        out <- lapply(out[[2]], function(x) {
            class(x) <- 'fedreg_document'
            return(x)
        })
    }
    return(out)
}

