get_meme <- 
function(id,
         ...){
    if(inherits(id, 'imgur_image') || inherits(id, 'imgur_gallery_image'))
        id <- id$id
    out <- imgurGET(paste0('gallery/g/memes/', id), ...)
    structure(out, class = 'imgur_gallery_image') # check this
}
