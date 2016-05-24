get_gallery <- 
function(section = 'hot',
         sort = 'viral',
         page = 0,
         window = 'day',
         showViral = TRUE,
         ...){
    stopifnot(section %in% c('hot', 'top', 'user'))
    if(!is.null(window))
       stopifnot(window %in% c('day', 'week', 'month', 'year', 'all'))
    stopifnot(sort %in% c('viral', 'time'))
    stopifnot(is.numeric(as.numeric(page)))
    stopifnot(is.logical(showViral))
    out <- imgurGET(paste0('gallery/', section, '/', sort,
                           ifelse(section == 'top', paste0('/', window), 
                                                    ''),
                           '/', page,
                           ifelse(section == 'user', paste0('?showViral=', showViral),
                                                     '')),
                    ...)
    out <- lapply(out, `class<-`, 'imgur_gallery_image')
    structure(out, class = 'imgur_gallery_album')
}
