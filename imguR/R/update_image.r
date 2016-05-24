update_image <- 
function(id, 
         title = NULL,
         description = NULL,
         ...){
    if(inherits(id, 'imgur_image'))
        id <- id$id
    out <- imgurPOST(paste0('image/', id, '/'), 
                     body = list(title = title, 
                                 description = description),
                     ...)
    structure(out, class = 'imgur_basic')
}
