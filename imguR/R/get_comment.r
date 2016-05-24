get_comment <-
function(comment,
         replies = FALSE,
         ...){
    if(replies)
        out <- imgurGET(paste0('comment/', comment, '/replies'), ...)
    else
        out <- imgurGET(paste0('comment/', comment), ...)
    structure(out, class = 'imgur_comment')
}
