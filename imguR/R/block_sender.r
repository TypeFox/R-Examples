block_sender <-
function(username,
         ...){
    if(!"token" %in% names(list(...)))
        stop("This operation can only be performed using an OAuth token.")
    out <- imgurPOST(paste0('conversations/block/', username), ...)
    structure(out, class = 'imgur_basic')
}
