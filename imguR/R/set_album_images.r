set_album_images <- 
function(album, 
         id,
         ...){
    if(!"token" %in% names(list(...)))
        stop("This operation can only be performed using an OAuth token.")
    if(inherits(album, 'imgur_album'))
        album <- album$id
    if(is.list(id))
        id <- sapply(id, `$`, 'id')
    b <- list(ids = paste(id, collapse = ','))
    out <- imgurPOST(paste0('album/', album), body = b, ...)
    structure(out, class = 'imgur_basic')
}