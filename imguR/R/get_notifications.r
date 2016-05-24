get_notifications <-
function(id = NULL,
         only_new = FALSE,
         ...){
    if(!"token" %in% names(list(...)))
        stop("This operation can only be performed using an OAuth token.")
    if(!is.null(id)) {
        out <- imgurGET(paste0('notification/', id), ...)
    } else {
        out <- imgurGET('notification/', 
                        body = list(new = only_new), ...)
    }
    structure(out, class = 'imgur_notification')
}
