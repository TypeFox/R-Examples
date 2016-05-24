report_comment <-
function(comment,
         ...){
    out <- imgurPOST(paste0('comment/', comment, '/report'), ...)
    structure(out, class = 'imgur_basic')
}
