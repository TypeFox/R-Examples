`combine` <-
function (pimgs, subset = 1:length(pimgs)) 
{
    bad <- unlist(lapply(pimgs, function(x) is.null(x$image)))[subset]
    subset <- subset[!bad]
    res <- as.image.pimg(pimgs[[subset[1]]])
    if (length(subset) == 1) 
        return(res)
    for (i in subset[-1]) {
        img <- pimgs[[i]]
        Xpos <- img$offset[1]
        Ypos <- img$offset[2]
        Xind <- Xpos:(Xpos + dim(img$image)[1] - 1)
        Yind <- Ypos:(Ypos + dim(img$image)[2] - 1)
        res$z[Xind, Yind] <- res$z[Xind, Yind] + img$image
    }
    res
}

