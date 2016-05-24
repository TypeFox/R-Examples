`DynamicList` <-
function(C,B,D) {
    L <- PathList(C,D) + 1
    bounds <- B[L]
    return(bounds)
}

