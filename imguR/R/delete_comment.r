delete_comment <-
function(comment, ...){
    if(!"token" %in% names(list(...)))
        stop("This operation can only be performed using an OAuth token.")
    out <- imgurDELETE(paste0('comment/', comment), ...)
    structure(out, class = 'imgur_basic')
}
