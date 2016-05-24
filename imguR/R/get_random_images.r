get_random_images <- 
function(page = 0, ...){
    stopifnot(is.numeric(as.numeric(page)))
    out <- imgurGET(paste0('gallery/random/random/',
                           ifelse(!is.null(page), page, NULL)),
                    ...)
    out <- lapply(out, `class<-`, 'imgur_gallery_image')
    structure(out, class = 'imgur_gallery_album')
}
