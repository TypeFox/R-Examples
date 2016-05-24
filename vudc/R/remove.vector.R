remove.vector <-
function (x, remove.absolute, remove.ratio) 
{
    x.removed <- x
    if (!is.na(remove.absolute)) {
        x.removed <- x[abs(x) < remove.absolute]
    }
    else if (!is.na(remove.ratio)) {
        x.lowerquantile <- quantile(x, remove.ratio)
        x.upperquantile <- quantile(x, 1 - remove.ratio)
        x.removed <- x[(x.lowerquantile[[1]] < x) & (x < x.upperquantile[[1]])]
    }
    return(x.removed)
}
