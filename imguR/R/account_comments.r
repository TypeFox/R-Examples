account_comments <-
function(account = 'me',
         ids = FALSE,
         ...){
    if(!"token" %in% names(list(...)) && account == 'me')
        stop("This operation can only be performed for account 'me' using an OAuth token.")
    if(ids) {
        out <- imgurGET(paste0('account/', account, '/comments/ids'), ...)
        structure(out, class = 'imgur_basic')
    } else {
        out <- imgurGET(paste0('account/', account, '/comments'), ...)
        structure(out, class = 'imgur_comment')
    }
}
