fr_agencies <- function(id=NULL, version='v1', ...){
    if(is.null(id))
        baseurl <- paste('https://www.federalregister.gov/api/',version,'/agencies', sep='')
    else
        baseurl <- paste('https://www.federalregister.gov/api/',version,'/agencies/', id, sep='')
    
    args <- NULL
    
    r <- GET(url = paste(baseurl, args, sep=''), ...)
    stop_for_status(r)
    response <- content(r, "text")
    out <- fromJSON(response)
    if(length(id)==1){
        class(out) <- 'fedreg_agency'
    } else if(length(out)>1){
        out <- lapply(out, function(x) {
            class(x) <- 'fedreg_agency'
            return(x)
        })
    }
    return(out)
}

print.fedreg_agency <- function(x, ...){
    if(!is.null(x$short_name))
        cat(x$short_name, ': ', sep='')
    cat(x$name, '\n', sep='')
    cat('ID:', x$id)
    if(!is.null(x$parent_id))
        cat(' Parent ID:',x$parent_id, '\n')
    else
        cat('\n')
    cat(strwrap(x$description),'\n')
    cat('URL: ',x$url,'\n')
    cat('URL for recent articles: ', x$recent_articles_url, '\n')
    cat('\n')
    invisible(x)
}
