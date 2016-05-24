fle.plot <-
function(lca.obj, at = 65, ...) {
    if (any(lca.obj$fit$y>10)) stop('Check fitted rates! - could be num. of deaths')
    dat <- demogdata(lca.obj[[4]], array(0, dim(lca.obj[[4]])), lca.obj$age,
                     lca.obj$year, 'mortality', lca.obj$label, names(lca.obj)[4])
    obs.le <- life.expectancy(dat, age=at)
    fcast <- forecast(lca.obj, shift=T, ...)
    # quick fix for latest versions of forecast  (NULL type!)
    fcast$type <- "mortality"
    dat$rate[[1]] <- lca.obj$fit$y
    if (any(dat$rate[[1]] < 0)) dat$rate[[1]] <- exp(dat$rate[[1]])
    fit.le <- life.expectancy(dat, age=at)
    fcast.le <- life.expectancy(fcast, age=at)
    fcast.leb <- cbind(fcast$year, life.expectancy(fcast, 'lower', age=at),
                       life.expectancy(fcast, 'upper', age=at))
    ind <- or(is.na(fcast.leb[,2]), is.na(fcast.leb[,3]))
    fcast.leb <- fcast.leb[!ind,]
    xlim <- range(dat$year, fcast$year, na.rm=T)
    ylim <- range(obs.le, fcast.leb[,2:3], na.rm=T)
    plot(obs.le, ylab=substitute(le[x], list(x=at)), xlab='Year', xlim=xlim, ylim=ylim)
    lines(fit.le, col='green', lty=3)
    filled.gap(fcast.leb, edges=F, append=TRUE)
    lines(fcast.le, col='red', lty=2)
    legend(coord('UL'), legend=c('Obs', 'Fit', 'Fcast'), lty=c(1, 3, 2),
           col=c(1, 'green', 'red'))
    title(paste('Forecasts of Life Expectancy at age', at))
    title(paste(lca.obj$lab, names(lca.obj)[4], sep=' : '), line=0.5,
          cex=0.6, font=4)
}
