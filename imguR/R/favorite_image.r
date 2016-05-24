favorite_image <- 
function(id, ...){
    if(inherits(id, 'imgur_image'))
        id <- id$id
    out <- imgurPOST(paste0('image/', id, '/favorite'), ...)
    structure(out, class = 'imgur_basic')
}