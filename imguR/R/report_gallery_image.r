report_gallery_image <-
function(id,
         ...){
    out <- imgurPOST(paste0('gallery/image/', id, '/report'), ...)
    structure(out, class = 'imgur_vote')
}
