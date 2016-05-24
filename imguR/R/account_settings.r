account_settings <-
function(bio = NULL,
         public_images = NULL,
         messaging_enabled = NULL,
         album_privacy = NULL,
         accepted_gallery_terms = NULL,
         ...){
    if(!"token" %in% names(list(...)))
        stop("This operation can only be performed using an OAuth token.")
    a <- list()
    if(!is.null(bio))
        a$bio <- bio
    if(!is.null(public_images))
        a$public_images <- public_images
    if(!is.null(messaging_enabled))
        a$messaging_enabled <- messaging_enabled
    if(!is.null(album_privacy)) {
        stopifnot(album_privacy %in% c('public', 'hidden', 'secret'))
        a$album_privacy <- album_privacy
    }
    if(!is.null(accepted_gallery_terms))
        a$accepted_gallery_terms <- accepted_gallery_terms
    if(any(!sapply(a, is.null))) {
        out <- imgurPOST(paste0('account/me/settings'), body = a, ...)
        structure(out, class = 'imgur_basic')
    } else {
        out <- imgurGET(paste0('account/me/settings'), ...)
        structure(out, class = 'imgur_account_settings')
    }
}
