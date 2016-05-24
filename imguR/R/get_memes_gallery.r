get_memes_gallery <- 
function(sort = 'viral',
         page = 0,
         window = NULL,
         ...){
    if(!is.null(window))
       stopifnot(window %in% c('day', 'week', 'month', 'year', 'all'))
    stopifnot(sort %in% c('viral', 'time', 'top'))
    stopifnot(!is.numeric(as.numeric(page)))
    out <- imgurGET(paste0('gallery/g/memes/', sort,
                           ifelse(!is.null(window), paste0('/', window, '/'), 
                                                    '/'),
                           page),
                    ...)
    out <- lapply(out, `class<-`, 'imgur_gallery_image')
    structure(out, class = 'imgur_gallery_album')
}
