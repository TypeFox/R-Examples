matfle.plot <-
function(lca.obj, lca.base, at = 65, label=NULL, ...) {
  age <- lca.obj[[1]]$age
  year <- lca.obj[[1]]$year
  dat <- lapply(lca.obj, function(x) demogdata(x[[4]], array(0, dim(x[[4]])),
          age, year, 'mortality', x$label, names(x)[4]))
  for(i in seq(dat)) {
    if (any(lca.obj[[i]]$fit$y>10))
        stop('Check fitted rates! - could be num. of deaths')
      else dat[[i]]$rate[[1]] <- lca.obj[[i]]$fit$y
    if (any(dat[[i]]$rate[[1]] < 0)) dat[[i]]$rate[[1]] <- exp(dat[[i]]$rate[[1]])
  }
  fit.le <- as.data.frame(lapply(dat, life.expectancy, age=at))
  fcast <- lapply(lca.obj, forecast, shift=T, ...)
  # quick fix for latest versions of forecast  (NULL type!)
  fcast <- lapply(fcast, function(x) {x$type <- "mortality"; x})
  fcast.le <- as.data.frame(lapply(fcast, life.expectancy, age=at))
  if (!missing(lca.base)) {
    if (any(age=!lca.base$age))
      stop('Mismatch of age range between base model and subsets!')
    if (any(lca.base$fit$y>10))
      stop('Check fitted rates of base model! - could be num. of deaths')
    # browser()
    dat <- demogdata(lca.base[[4]], array(0, dim(lca.base[[4]])), lca.base$age,
          lca.base$year, 'mortality', lca.base$label, names(lca.base)[4])
    dat$rate[[1]] <- lca.base$fit$y
    if (any(dat$rate[[1]] < 0)) dat$rate[[1]] <- exp(dat$rate[[1]])
    fit.le$base <- life.expectancy(dat, age=at)
    fit.le <- fit.le[c(ncol(fit.le), 1:(ncol(fit.le)-1))]
    fcast <- forecast(lca.base, shift=T, ...)
    # quick fix for latest versions of forecast  (NULL type!)
    fcast$type <- "mortality"
    fcast.le$base <- life.expectancy(fcast, age=at)
    fcast.le <- fcast.le[c(ncol(fcast.le), 1:(ncol(fcast.le)-1))]
  }
  fyear <- as.numeric(dimnames(fcast.le)[[1]])
  xlim <- xaxp <- range(year, fyear)
  xlim[2] <- round(xlim[2]+0.1*diff(xlim))
  xlim[1] <- round(xlim[1]-0.35*diff(xlim))
  ylim <- range(fit.le, fcast.le, na.rm=T)
  # get the xaxis marks (a bit overcomplicated, but it looks good)
  xaxp <- round(xaxp)+5
  xaxp <- xaxp-mod(xaxp, 5)
  xaxp <- c(xaxp, length(seq(xaxp[1], xaxp[2], by=5))-1)
  matplot(year, fit.le, ylab=substitute(le[x], list(x=at)), xlab='Year',
    lty=seq(fit.le), xlim=xlim, ylim=ylim, type='l', col=rainbow(ncol(fit.le)),
    xaxp=xaxp)
  points(year, fit.le$base, pch='+', cex=0.75, col='red')
  matlines(fyear, fcast.le, col=rainbow(ncol(fcast.le)), lty=seq(fit.le), )
  points(fyear, fcast.le$base, pch='+', cex=0.75, col='red')
  text(xlim[2]-0.05*diff(xlim), fcast.le[nrow(fcast.le),], label=names(fcast.le),
     col=rainbow(ncol(fcast.le)), adj=0, cex=0.7)
  legend(coord(c(0.015,0.7)), legend=names(fit.le), lty=seq(fit.le), col=rainbow(ncol(fit.le)),
    text.col=rainbow(ncol(fit.le)), pch=c('+', rep('', length(lca.obj))), y.intersp=0.95,
    cex=0.95, yjust=0.5)
  title(paste('Forecasts of Life Expectancy at age', at))
  title(label, line=0.5, cex=0.6, font=4)
}
