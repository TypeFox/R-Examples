get_reddit_gallery <- 
function(subreddit,
         sort = 'time',
         page = 0,
         window = 'day',
         ...){
    if(!is.null(window))
       stopifnot(window %in% c('day', 'week', 'month', 'year', 'all'))
    stopifnot(sort %in% c('time', 'top'))
    stopifnot(is.numeric(as.numeric(page)))
    out <- imgurGET(paste0('gallery/r/', 
                           paste0(subreddit, '/'),
                           sort, '/',
                           ifelse(sort == 'top', paste0(window), 
                                                 ''),
                           '/', page),
                    ...)
    out <- lapply(out, `class<-`, 'imgur_gallery_image')
    structure(out, class = 'imgur_gallery_album')
}
