print.rhdata <-
function(x, ...){
    cat('Multidimensional Mortality data for:', x$label, sqb(x$name), '\n')
    cat('Across covariates:')
    cat(paste("\n\t years:", min(x$year), "-", max(x$year)))
    cat(paste("\n\t ages: ", min(x$age), "-", max(x$age)))
    for(i in seq(length(x$covariates))){
        cat(paste('\n\t ', names(x$covariates)[i], ': ',
                  paste(x$covariates[[i]], collapse=', '), sep=''))
    }
    cat('\n')   
}
