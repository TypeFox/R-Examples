delete_conversation <-
function(id, ...){
    out <- imgurDELETE(paste0('conversations/', id), ...)
    structure(out, class = 'imgur_basic')
}
