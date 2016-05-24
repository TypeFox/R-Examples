summary.profLinear <-
function(object, ...) {
    r <- length(unique(object$clust))
    q <- length(object$m[[1]])
    uclust <- object$clust[match(unique(object$group), object$group)]
    tclust <- table(uclust)
    tclusto <- table(object$clust)
    term.labels <- attr(terms(object$model), 'term.labels')
    if(attr(terms(object$model), 'intercept'))
        term.labels <- c('(intercept)', term.labels)
    cat('----------\n')
    out <- vector('list', length=r)
    for(cls in 1:r) {
        out[[cls]]$groups <- tclust[cls]
        out[[cls]]$obs <- tclusto[cls]
        cat('cluster: ', cls, '\n', sep='')        
        cat('groups:  ', tclust[cls], '\n', sep='') 
        cat('obs:     ', tclusto[cls], '\n', sep='')
        #Laplace approximation to 95%CI for coefficients
        M <- object$m[[cls]]
        S <- 1/sqrt(diag(object$s[[cls]] * (object$a[[cls]]+q/2) /
             object$b[[cls]])) 
        lower <- M - 1.92 * S
        upper <- M + 1.92 * S
        info <- data.frame(estimate=object$m[[cls]], 
            lower95=lower, upper95=upper)
        row.names(info) <- term.labels
        out[[cls]]$summary <- info
        print(format(info, nsmall=3, digits=3))
        cat('----------\n')
    }
    invisible(out)
}

