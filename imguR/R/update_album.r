update_album <- 
function(album, 
         id = NULL,
         title = NULL,
         description = NULL,
         privacy = NULL,
         layout = NULL,
         cover_id = NULL,
         ...){
    if(inherits(album, 'imgur_album'))
        album <- album$id
    b <- list()
    if(!is.null(title))
        b$title <- title
    if(!is.null(description))
        b$description <- description
    if(!is.null(privacy)) {
        stopifnot(privacy %in% c('public', 'hidden', 'secret'))
        b$privacy <- privacy
    }
    if(!is.null(layout)) {
        stopifnot(layout %in% c('blog', 'grid', 'horizontal', 'vertical'))
        b$layout <- layout
    }
    if(!is.null(cover_id)){
        if(inherits(cover_id, 'imgur_image'))
            cover_id <- cover_id$id
        b$cover_id <- cover_id
    }
    if(!is.null(id)){
        if(is.list(id))
            id <- sapply(id, `$`, 'id')
        b$ids <- paste(id, collapse = ',')
    }
    out <- imgurPOST(paste0('album/', album), body = b, ...)
    structure(out, class = 'imgur_basic')
}
