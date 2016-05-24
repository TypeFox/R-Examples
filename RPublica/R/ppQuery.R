ppQuery <-
function(op, baseurl, args=NULL, ...){
    response <- GET(paste(baseurl, op, '.json', args, sep=''), ...)
    stop_for_status(response)
    fromJSON(content(response, as = 'text'), flatten = TRUE)
}
