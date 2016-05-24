"getstnname" <-
function(id) {
    # Lookup MSC station data from internal data.frame
    name <- seas::mscstn$name[match(id, seas::mscstn$nid)]
    if (length(name) == 0 || is.na(name))
        name <- NULL
    name
}
