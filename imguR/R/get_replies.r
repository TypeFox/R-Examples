get_replies <-
function(only_new = FALSE, ...){
    if(!"token" %in% names(list(...)))
        stop("This operation can only be performed using an OAuth token.")
    out <- imgurGET(paste0('account/me/notifications/replies'),
                    body = list(new = only_new), ...)
    structure(out, class = 'imgur_notification')
}
