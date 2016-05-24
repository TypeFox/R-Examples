add_gallery_images <- 
function(album = NULL, 
         id = NULL,
         title = NULL,
         ...){
    b <- list(title = title,
              terms = 1)
    if(!is.null(image)) {
        if(is.null(title))
            b$title <- get_image(id, ...)$title
        out <- imgurPOST(paste0('gallery/image/', image), body = b, ...)
        structure(out, class = 'imgur_gallery_image')
    }
    if(!is.null(album)){
        if(is.null(title))
            b$title <- get_album(album, ...)$title
        out <- imgurPOST(paste0('gallery/album/', album), body = b, ...)
        structure(out, class = 'imgur_gallery_album')
    }
}