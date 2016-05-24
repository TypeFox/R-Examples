geos <- function(state=NULL, ...){
    if(is.null(state))
        out <- ppQuery('geos', baseurl='http://projects.propublica.org/forensics/', ...)
    else {
        out <- ppQuery(paste('geos', state, sep='/'),
                       baseurl='http://projects.propublica.org/forensics/', ...)
    }
    return(out)
}

systems <- function(id, ...){
    out <- ppQuery(paste('systems', id, sep='/'),
                   baseurl='http://projects.propublica.org/forensics/', ...)
    return(out[[1]])
}
