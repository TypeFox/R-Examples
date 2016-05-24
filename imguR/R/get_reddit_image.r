get_reddit_image <- 
function(subreddit,
         id,
         ...){
    if(inherits(id, 'imgur_image') || inherits(id, 'imgur_gallery_image'))
        id <- id$id
    out <- imgurGET(paste0('gallery/r/', subreddit, '/', id), ...)
    structure(out, class = 'imgur_gallery_image')
}
