get_album <- 
function(album, ...){
    if(inherits(album, 'imgur_album') || inherits(album, 'imgur_gallery_album'))
        album <- album$id
    out <- imgurGET(paste0('album/', album), ...)
    out[['images']] <- lapply(out[['images']], `class<-`, 'imgur_image')
    structure(out, class = 'imgur_album')
}