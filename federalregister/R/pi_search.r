pi_search <- function(..., fields=NULL, per_page=NULL, page=NULL, version='v1', getopts = NULL) {
    
    baseurl <- paste('https://www.federalregister.gov/api/',version,
                     '/public-inspection-documents.json?', sep='')

    query <- list(...)
    
    # api docs mention separate by-date API, but same as search:
    # conditions[available_on]=date
    if('date' %in% names(query))
        names(query)[names(query)=='date'] <- 'available_on'
    query <- curl_escape(paste('conditions[',names(query),']=',query,sep='',collapse='&'))
    
    # handle pagination
    if(!is.null(per_page) && as.numeric(per_page)>1000)
        stop("'per_page' cannot be greater than 1000")
    else if(!is.null(per_page) & !is.null(page))
        p <- paste('per_page=',per_page,'&page=',page,sep='')
    else if(!is.null(per_page) & is.null(page))
        p <- paste('per_page=',per_page,sep='')
    else if(!is.null(page))
        p <- paste('page=',page,sep='')
    else 
        p <- NULL
    
    if(!is.null(fields)){
        fields <- paste(curl_escape(paste('fields[]',fields,sep='=')), collapse='&')
        args <- paste(fields,query,p,sep='&')
    } else
        args <- paste(query,p,sep='&')
    
    r <- do.call("GET", c(list(url = paste(baseurl, args, sep='')), getopts))
    stop_for_status(r)
    response <- content(r, "text")
    out <- fromJSON(response)
    out$results <- lapply(out$results, `class<-`, 'fedreg_document')
    return(out)
}
