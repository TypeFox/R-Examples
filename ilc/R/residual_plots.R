residual_plots <-
function(lca.obj, file=paste('res', deparse(substitute(lca.obj)), 
                                               'ps', sep='.'), view=T){
    if (!is.null(file)) postscript(file, horizontal=FALSE)
    else old.par <- par(no.readonly=T)
    par(mfrow=c(3,1), mar=c(5, 5, 4, 2) + 0.1)
    ylab <- paste(lca.obj$res$type, 'residuals')
    # a) by calendar year:
    matplot(lca.obj$year, t(lca.obj$res$y), xlab='calendar year', ylab=ylab, 
            pch=3, cex=0.5, col=1)
    title(paste(lca.obj$res$ynam, 'by year'))
    title(paste(lca.obj$lab, names(lca.obj)[4], sep=' : '), line=0.5, cex.main=0.8)
    # b) by age
    matplot(lca.obj$age, lca.obj$res$y, xlab='age', ylab=ylab, pch=3, cex=0.5, col=1)
    title(paste(lca.obj$res$ynam, 'by age'))
    title(paste(lca.obj$lab, names(lca.obj)[4], sep=' : '), line=0.5, cex.main=0.8)
    # c) by cohorts
    coh <- cohort(lca.obj$year, lca.obj$age)
    # plot(coh, lca.obj$res$y, type='n', xlim=c(1885,1945), )
    plot(coh, lca.obj$res$y,  , xlab='year of birth', ylab=ylab, pch=3, cex=0.5, 
         col=1)
    title(paste(lca.obj$res$ynam, 'by cohort'))
    title(paste(lca.obj$lab, names(lca.obj)[4], sep=' : '), line=0.5, cex.main=0.8)
    if (!is.null(file)) { dev.off(); if (view) ps.view(file) }
    else invisible(par(old.par))
}
