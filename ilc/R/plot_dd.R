plot_dd <-
function(dd.obj, year=dd.obj$year, col=rainbow(length(year), start=0.1),
      lpos='UL', lpar=list(), ppar = T, ...){
  tmp <- year%in%dd.obj$year
  if (!all(tmp)) year <- year[tmp]
  if (dd.obj$type == 'deaths') {
    dd.obj$type <- 'mortality'
  }
  plot(dd.obj, year=year, col=col, ...)

  def.par <- list(x=coord(lpos))
  tmp <- pmatch(names(lpar), 'title', nomatch=F)
  if (!any(tmp)) def.par$title <- 'Year'
  lpar <- c(lpar, def.par)
  lpar$text.col <- col
  lpar$col <- col
  lpar$legend <- year
  if (ppar){
      dots <- list(...)
      tmp <- c('lty', 'lwd', 'pch')
      tmp <- names(dots)%in%tmp
      if (any(tmp)) lpar <- c(lpar, dots[tmp])
  }
  do.call('legend', lpar)
}
