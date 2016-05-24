add_album_images <- 
function(album, 
         id,
         ...){
    if(!"token" %in% names(list(...)))
        stop("This operation can only be performed using an OAuth token.")
    if(inherits(album, 'imgur_album'))
        album <- album$id
    if(inherits(id, 'imgur_image')) {
        id <- id$id
    } else if(is.list(id)) {
        id <- sapply(id, `$`, 'id')
    }
    b <- list(ids = paste(id, collapse = ','))
    out <- imgurPOST(paste0('album/', album, '/add'), body = b, ...)
    structure(out, class = 'imgur_basic')
}