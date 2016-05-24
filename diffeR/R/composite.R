composite <- function (comp, ref, factor){

    nRefCats <- length(unique(ref))
    nCompCats <- length(unique(comp))

    if (nCompCats != nRefCats) 
        stop("number of layers in input rasters must be equal")

    conf.matrix <- matrix(nrow = nCompCats, ncol = nRefCats)

    compAg <- memberships(comp, fact=factor)
    refAg <- memberships(ref, fact=factor)

    weights <- comp-comp+1

#    comp.lys <- layerize(comp)
#    compAg <- aggregate(comp.lys, fact=factor, fun='mean')
    
#    ref.lys <- layerize(ref)
#    refAg <- aggregate(ref.lys, fact=factor, fun='mean')

    wghtAg <- aggregate(weights, fact=factor, fun=sum)

    sumDiagPix <- 0
    for (i in 1:nlayers(refAg)) {
        min.ii <- min(compAg[[i]], refAg[[i]])
        sumDiagPix <- sumDiagPix + min.ii
    }
    for (i in 1:nlayers(compAg)) {
        minsi <- min(compAg[[i]], refAg[[i]])
        for (j in 1:nlayers(refAg)) {
            minsj <- min(compAg[[j]], refAg[[j]])
            Aij <- wghtAg*((compAg[[i]] - minsi) * (refAg[[j]] - minsj))/(1 - 
                sumDiagPix)
            sum.Aij <- cellStats(Aij, stat = "sum")
            conf.matrix[i, j] <- sum.Aij/cellStats(wghtAg, stat = "sum")
        }
    }
    for (i in 1:nlayers(refAg)) {
        min.ii <- min(compAg[[i]], refAg[[i]])*wghtAg
        sum.min.ii <- cellStats(min.ii, stat = "sum")
        conf.matrix[i, i] <- sum.min.ii/cellStats(wghtAg, stat = "sum")
    }
    colnames(conf.matrix) <- paste("ref.", unique(ref), sep="")
    rownames(conf.matrix) <- paste("comp.", unique(comp), sep="")
    return(conf.matrix)
}