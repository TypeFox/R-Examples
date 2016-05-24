T.maxmean <-
function(geneZ){

    Zvec = geneZ[,2]
    gname0 = as.vector(geneZ[,1])
    names(Zvec) = gname0

    Tvec <- tapply(Zvec, gname0, maxmean.op)
    mvec <- tapply(gname0, gname0, length)
    Txx = cbind(Tvec,mvec)
    colnames(Txx) = c('T','m')
    rownames(Txx) = names(mvec)
    return(Txx)
    }
