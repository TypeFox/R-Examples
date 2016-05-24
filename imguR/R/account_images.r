account_images <-
function(account = 'me',
         page = NULL, 
         ids = TRUE,
         ...){
    if(!"token" %in% names(list(...)) && account == 'me')
        stop("This operation can only be performed for account 'me' using an OAuth token.")
    if(ids) {
        out <- imgurGET(paste0('account/', account, '/images/ids'), ...)
        structure(out, class = 'imgur_basic')
    } else {
        out <- imgurGET(paste0('account/', account, '/images/', page), ...)
        structure(out, class = 'imgur_image') # check this
    }
    
}
