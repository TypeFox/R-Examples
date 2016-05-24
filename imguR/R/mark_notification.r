mark_notification <-
function(id,
         ...){
    if(!"token" %in% names(list(...)))
        stop("This operation can only be performed using an OAuth token.")
    out <- imgurPOST(paste0('notification/', id), ...)
    structure(out, class = 'imgur_basic')
}
