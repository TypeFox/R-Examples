search_gallery <- 
function(query,
         sort = 'time',
         page = 1,
         ...){
    stopifnot(sort %in% c('time', 'top'))
    stopifnot(is.numeric(as.numeric(page)))
    out <- imgurGET(paste0('gallery/search/', sort,
                           ifelse(!is.null(page), paste0('/', page, '/'), 
                                                    '/'),
                           '?q=', query),
                    ...)
    out <- lapply(out, `class<-`, 'imgur_gallery_image')
    structure(out, class = 'imgur_gallery_album')
}
