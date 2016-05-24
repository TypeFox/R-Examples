rmProfiles.wprof <-
function(y, v, ...) {
    if (is.logical(v)) v <- which(v)
    if (is.character(v)) v <- which(rownames(prof$profiles) %in% v)
    prof$profiles <- prof$profiles[-v,]
    prof$freq <- prof$freq[-v]
    return(prof)
}
