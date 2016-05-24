flc.plot <-
function(lca.obj, at = 65, ...) {
    p.old <- par()
    par(mfrow=c(1,2), mar=c(5, 5, 4, 2) + 0.1)
    fcast <- forecast(lca.obj, , shift=F, ...)    # for older forecast versions: 
    # quick fix for latest versions of forecast  (NULL type!)
    fcast$type <- "mortality"
    plot(fcast$kt, xlab='Year', ylab=substitute(kappa[t][' '] (a), list(a=lca.obj$adj)))
    title(paste(lca.obj$lab, names(lca.obj)[4], sep=' : '), line=0.5,
          cex=0.6, font=4)
    fle.plot(lca.obj, at=at, ...)
    invisible(par(p.old))
}
