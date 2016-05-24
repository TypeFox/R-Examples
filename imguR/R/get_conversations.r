get_conversations <-
function(id = NULL, ...){
    if(!is.null(id)) {
        out <- imgurGET(paste0('conversations/', id), ...)
    } else {
        out <- imgurGET('conversations/', ...)
    }
    structure(out, class = 'imgur_message')
}
