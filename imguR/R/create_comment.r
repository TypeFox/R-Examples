create_comment <-
function(id,
         comment,
         parent = NULL,
         ...){
    if(!"token" %in% names(list(...)))
        stop("This operation can only be performed using an OAuth token.")
    if(inherits(id, 'imgur_image'))
        id <- id$id
    b <- list(image_id = id,
              comment = comment)
    if(!is.null(parent))
        b$parent_comment <- parent
    out <- imgurPOST('comment/',
                     body = b,
                     ...)
    structure(out, class = 'imgur_basic')
}
