account_favorites <-
function(account = 'me',
         gallery = FALSE, 
         ...){
    if(!"token" %in% names(list(...)) && account == 'me')
        stop("This operation can only be performed for account 'me' using an OAuth token.")
    out <- imgurGET(paste0('account/', account, 
                           ifelse(gallery, '/gallery_favorites', '/favorites')),
                    ...)
    if(gallery)
        structure(out, class = 'imgur_gallery_album') # check this
    else
        structure(out, class = 'imgur_image') # check this
}
