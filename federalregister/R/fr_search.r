fr_search <- function(..., fields=NULL, per_page=NULL, page=NULL,
                      order='relevance', version='v1', getopts = NULL) {

    baseurl <- paste('https://www.federalregister.gov/api/',version,
                     '/articles.json?', sep='')
    
    query <- list(...)
    
    # need to process `publication_date` list
    if('publication_date' %in% names(query)){
        w <- which(names(query)=='publication_date')
        p <- query$publication_date
        names(p) <- paste('publication_date][',names(p),sep='')
        query <- query[-w]
        query <- c(query, p)
    }
    # need to process `effective_date` list
    if('effective_date' %in% names(query)){
        w <- which(names(query)=='effective_date')
        p <- query$effective_date
        names(p) <- paste('effective_date][',names(p),sep='')
        query <- query[-w]
        query <- c(query, p)
    }
    # need to process `cfr` list
    if('cfr' %in% names(query)){
        w <- which(names(query)=='cfr')
        p <- query$cfr
        names(p) <- paste('cfr][',names(p),sep='')
        query <- query[-w]
        query <- c(query, p)
    }
    # need to process `near` list
    if('near' %in% names(query)){
        w <- which(names(query)=='near')
        p <- query$near
        names(p) <- paste('near][',names(p),sep='')
        query <- query[-w]
        query <- c(query, p)
    }
    
    query <- paste(curl_escape(paste('conditions[',names(query),']=',query,sep='')),collapse='&')
    
    # handle pagination
    if(!is.null(per_page) && as.numeric(per_page)>1000) {
        stop("'per_page' cannot be greater than 1000")
    } else if(!is.null(per_page) & !is.null(page)) {
        p <- paste('per_page=',per_page,'&page=',page,sep='')
    } else if(!is.null(per_page) & is.null(page)) {
        p <- paste('per_page=',per_page,sep='')
    } else if(!is.null(page)) {
        p <- paste('page=',page,sep='')
    } else {
        p <- NULL
    }
    
    if(!is.null(fields)){
        fields <- paste(paste(curl_escape('fields[]'),fields,sep='='), collapse='&')
        args <- paste(fields,query,p,sep='&')
    } else {
        args <- paste(query,p,sep='&')
    }
    
    
    r <- do.call("GET", c(list(url = paste(baseurl, args, sep='')), getopts))
    stop_for_status(r)
    response <- content(r, "text")
    out <- fromJSON(response)
    out$results <- lapply(out$results, `class<-`, 'fedreg_document')
    return(out)
}
