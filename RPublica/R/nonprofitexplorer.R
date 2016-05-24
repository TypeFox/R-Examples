np_search <- function(q=NULL, order=NULL, sort=NULL, state=NULL, ntee=NULL, subsection=NULL, ...){
    args <- formals()
    if(!is.null(q))
        args[names(args)=='q'] <- curl_escape(q)
    if(!is.null(order))
        args[names(args)=='q'] <- curl_escape(order)
    names(args)[names(args)=='sort'] <- curl_escape('sort_order')
    names(args)[names(args)=='state'] <- curl_escape('state[id]')
    names(args)[names(args)=='ntee'] <- curl_escape('ntee[id]')
    names(args)[names(args)=='subsection'] <- curl_escape('c_code[id]')
    args <- paste('?', paste(names(args), args, sep='=', collapse='&'), sep='')
    out <- ppQuery('search',
                   'https://projects.propublica.org/nonprofits/api/v1/',
                   args=args, ...)
    return(out)
}

np_org <- function(ein, ...){
    out <- ppQuery(paste('organizations',ein,sep='/'),
                   'https://projects.propublica.org/nonprofits/api/v1/', ...)
    return(out)
}
