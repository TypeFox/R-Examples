delete_image <- 
function(id, ...){
    if(inherits(id, 'imgur_image'))
        id <- id$id
    out <- imgurDELETE(paste0('image/', id), ...)
    structure(out, class = 'imgur_basic')
}