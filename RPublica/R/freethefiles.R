market <- function(slug=NULL, ...){
    if(is.null(slug))
        out <- ppQuery('markets', baseurl='https://projects.propublica.org/free-the-files/', ...)
    else
        out <- ppQuery(paste('markets', slug, sep='/'),
                       baseurl='https://projects.propublica.org/free-the-files/', ...)
    return(out)
}

station <- function(callsign, ...){
    out <- ppQuery(paste('stations', callsign, sep='/'),
                   baseurl='https://projects.propublica.org/free-the-files/', ...)
    return(out)
}

committee <- function(slug=NULL, ...){
    if(is.null(slug))
        out <- ppQuery('committees', 
                       baseurl='https://projects.propublica.org/free-the-files/', ...)
    else
        out <- ppQuery(paste('committees', slug, sep='/'),
                       baseurl='https://projects.propublica.org/free-the-files/', ...)
    return(out)
}

filing <- function(id, ...){
    out <- ppQuery(paste('filings', id, sep='/'), 
                   baseurl='https://projects.propublica.org/free-the-files/', ...)
    return(out[[1]])
}
