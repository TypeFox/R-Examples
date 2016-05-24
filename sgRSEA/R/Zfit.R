Zfit <-
function(dat, p.null){
    Zvec <- apply(dat[,3:4], MARGIN=1, Zstat, p.null=p.null)
    gname0 = as.vector(dat[,2])
    names(Zvec) = gname0
    return(Zvec)
    }
