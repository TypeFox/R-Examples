remove_gallery_image <- 
function(id,
         ...){
    if(!"token" %in% names(list(...)))
        stop("This operation can only be performed using an OAuth token.")
    out <- imgurDELETE(paste0('gallery/', id), ...)
    structure(out, class = 'imgur_basic')
}