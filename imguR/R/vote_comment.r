vote_comment <-
function(comment,
         vote = 'up',
         ...){
    if(!"token" %in% names(list(...)))
        stop("This operation can only be performed using an OAuth token.")
    if(!vote %in% c('up', 'down'))
        stop("'vote' can only be 'up' or 'down'")
    out <- imgurPOST(paste0('comment/', comment, '/vote/', vote), ...)
    structure(out, class = 'imgur_basic')
}
