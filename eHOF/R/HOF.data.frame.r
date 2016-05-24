"HOF.data.frame" <- function (
		veg, 
		grad, 
		M = max(veg), 
		freq.limit = 10, 
		...)
{
    if(names(veg)[1] == 'RELEVE_NR' & names(veg)[2] == 'SPECIES_NR') print('tvabund format detected, use tv.veg first')   else {
    if(!is.vector(grad)) print('Gradient must be a vector')
    freq <- apply(veg > 0, 2, sum)
    veg <- veg[, freq >= freq.limit, drop = FALSE]
    out <- lapply(veg, HOF.default, grad = grad, M = M, ...)
    for(i in seq_along(out)) out[[i]]$y.name <- names(veg)[i]
    class(out) <- "HOF.list"
    out
    }
}


# '[.HOF.list' <- function(x, ...) {
#    x <- UseMethod('[')
#    class(x) <- 'HOF.list'
# }
