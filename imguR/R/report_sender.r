report_sender <-
function(username,
         ...){
    out <- imgurPOST(paste0('conversations/report/', username), ...)
    structure(out, class = 'imgur_basic')
}
