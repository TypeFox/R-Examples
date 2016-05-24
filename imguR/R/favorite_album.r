favorite_album <- 
function(album, ...){
    if(inherits(album, 'imgur_album'))
        album <- album$id
    out <- imgurPOST(paste0('album/', album, 'favorite'), ...)
    structure(out, class = 'imgur_basic')
}