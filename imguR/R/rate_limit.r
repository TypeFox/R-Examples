rate_limit <-
function(...){
    out <- imgurGET('credits/', ...)
    out$UserReset <- as.POSIXct(as.numeric(out$UserReset), 
                                origin = '1970-01-01')
    structure(out, class = 'imgur_basic')
}
