matflc.plot <-
function(lca.obj, lca.base, at = 65, label=NULL, ...) {
  p.old <- par(no.readonly=T)
  par(mfrow=c(1,2), mar=c(5, 5, 4, 2) + 0.1)
  age <- lca.obj[[1]]$age
  year <- lca.obj[[1]]$year
  kt <- as.data.frame(lapply(lca.obj, function(x) x$kt))
  fcast <- as.data.frame(lapply(lca.obj,
                          function(x) forecast(x, shift=F, ...)$kt[[1]]))
  if (!missing(lca.base)){
    if (any(age!=lca.base$age))
      stop('Mismatch of age range between base model and subsets!')
    kt$base <- lca.base$kt
    kt <- kt[c(ncol(kt), 1:(ncol(kt)-1))]
    fcast$base <- forecast(lca.base, shift=F, ...)$kt[[1]]
    fcast <- fcast[c(ncol(fcast), 1:(ncol(fcast)-1))]
  }
  fyear <- year[length(year)]+seq(nrow(fcast))
  xlim <- range(year,  fyear)
  xlim[2] <- round(xlim[2]+0.1*diff(xlim))
  ylim <- range(kt, fcast, na.rm=T)
  matplot(year, kt, xlab='Year', ylab=substitute(kappa[t][' '] (a), list(a=lca.obj[[1]]$adj)), xlim=xlim,  ylim=ylim,
    type='l', lty=seq(fcast), col=rainbow(ncol(fcast)))
  matlines(fyear, fcast, lty=seq(fcast), col=rainbow(ncol(fcast)))
  if (!missing(lca.base)){
    points(year, kt$base, pch='+', cex=0.75, col='red')
    points(fyear, fcast$base, pch='+', cex=0.75, col='red')
  }
  # text(xlim[2]-0.07*diff(xlim), fcast[nrow(fcast),], label=names(fcast),
  #   col=rainbow(ncol(fcast)), adj=0, cex=0.8)
  legend(coord('UR'), legend=names(fcast), lty=seq(fcast), col=rainbow(ncol(fcast)),
    text.col=rainbow(ncol(fcast)), pch=c('+', rep('', length(lca.obj))), y.intersp=0.95)
  title('Forecasts of kt from Random walk with drift')
  title(label, line=0.5, cex=0.6, font=4)
  matfle.plot(lca.obj, lca.base, at=at, label=label, ...)
  invisible(par(p.old))
}
