fitted_plots <-
function(lca.obj, file=paste('fit', deparse(substitute(lca.obj)),
             'ps', sep='.'), view=T, labs=T, col){
  if (!is.null(file)) postscript(file, horizontal=T)
    else old.par <- par(no.readonly=T)
  par(mfrow=c(1,2), mar=c(5, 5, 4, 2) + 0.1)
  ylab <- paste('fitted', lca.obj$fit$type)
  nyears <- length(lca.obj$year)
  age <- lca.obj$age
  year <- lca.obj$year
  # a) by calendar year:
  xlim <- range(year)
  if (labs) xlim[2] <- xlim[2]+0.15*diff(xlim)
  if (missing(col)) col <- rainbow(nyears)
  matplot(year, t(lca.obj$fit$y), type='l', col=col,
    xlab='calendar year', ylab=ylab, xlim=xlim)
  if (labs){
    intv <- diff(age)
    if (any(intv > 1)) {
      warning('Applying age-intervals!')
      intv <- c(intv,1)
      intv <- apply(cbind(age, age+intv-1), 1, str.range)
      intv[length(intv)] <- paste(age[length(age)], '>', sep='')
      age <- intv
    }
    text(xlim[2]-1, lca.obj$fit$y[,nyears], labels=age, col=col, cex=0.9)
  }
  title(paste(lca.obj$fit$ynam, 'by year'))
  title(paste(lca.obj$lab, names(lca.obj)[4], sep=' : '), line=0.5, cex.main=0.8)
  # b) by age
  xlim <- range(age)
  if (labs) xlim[1] <- xlim[1]-0.25*diff(xlim)
  plot(lca.obj$fit, col=col, ylab=ylab, xlim=xlim)
  if (labs) legend(coord('UL'), legend=year, text.col=col,
      title='Years', x.intersp=-0.2, y.intersp=0.9, cex=0.9)
  title(paste(lca.obj$fit$ynam, 'by age'))
  title(paste(lca.obj$lab, names(lca.obj)[4], sep=' : '), line=0.5, cex.main=0.8)
  if (!is.null(file)) { dev.off(); if (view) ps.view(file) }
    else invisible(par(old.par))
}
