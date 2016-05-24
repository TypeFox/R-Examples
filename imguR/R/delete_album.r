delete_album <- 
function(album, ...){
    if(inherits(album, 'imgur_album'))
        album <- album$id
    out <- imgurDELETE(paste0('album/', album), ...)
    structure(out, class = 'imgur_basic')
}