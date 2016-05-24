get_album_images <- 
function(album, ...){
    if(inherits(album, 'imgur_album') || inherits(album, 'imgur_gallery_album'))
        album <- album$id
    out <- imgurGET(paste0('album/', album, '/images'), ...)
    structure(out, class = 'imgur_image')
}
