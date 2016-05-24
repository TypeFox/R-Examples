remove.list <-
function (x, remove.absolute, remove.ratio) 
{
    x.removed <- x
    if (!is.na(remove.absolute)) {
        x.removed <- list()
        for (i in 1:length(x)) {
            x.actual <- x[[i]]
            x.removed[[i]] <- x.actual[abs(x.actual) < remove.absolute]
        }
    }
    else if (!is.na(remove.ratio)) {
        x.composite <- c()
        for (i in 1:length(x)) {
            x.composite <- c(x.composite, x[[i]])
        }
        x.lowerquantile <- quantile(x.composite, remove.ratio)[[1]]
        x.upperquantile <- quantile(x.composite, 1 - remove.ratio)[[1]]
        x.removed <- list()
        for (i in 1:length(x)) {
            x.column <- x[[i]]
            x.removed[[i]] <- x.column[(x.lowerquantile < x.column) & 
                (x.column < x.upperquantile)]
        }
    }
    return(x.removed)
}
