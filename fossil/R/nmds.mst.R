`nmds.mst` <-
function(nmds,mst,...) {
    a <- nmds$stress
    b <- which.min(a)
    plot(nmds$conf[[b]][,1],nmds$conf[[b]][,2], ...)
    for (i in 1:ncol(mst)) {
      for (j in 1:ncol(mst)) {
        if (mst[i,j]==1) {lines(c(nmds$conf[[b]][i,1],nmds$conf[[b]][j,1]),c(nmds$conf[[b]][i,2],nmds$conf[[b]][j,2]),)}}}

}

