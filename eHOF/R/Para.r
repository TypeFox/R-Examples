"Para" <-  function(resp, ...) {
    UseMethod("Para")
}

"Para.HOF" <- function (
  resp, 
  model, 
  newdata = NULL, 
  ...)
{
  if (missing(model)) model <- pick.model(resp, gam=FALSE, ...)
  if (is.null(newdata)) x <- seq(-1, 2, length.out=10000) else x <- scale01(newdata, ...)
  M <- resp$M
  opt <- Para_opt(resp, model=model, newdata=x, ...)
  border <- Para_niche(resp, newdata=x, model=model, top=opt$top, optima=opt$opt, mini=opt$mini, pess=opt$pess, ...)
  slope <- Para_deriv(resp, newdata=x, p=resp$models[[model]]$par, model)
  infl <- Para_deriv(resp, newdata=x, p=resp$models[[model]]$par, model, optima=opt$opt, pessima=opt$pess, type='inflection')
  max.sl <- max(abs(slope))
  out <- list(species = resp$y.name, abund.sum = sum(resp$y/M), range = resp$range, model = model, para = resp$models[[model]]$par, M = M, mini = opt$mini, pess= opt$pess, top = opt$top, opt = opt$opt, max.slope=max.sl,  inflection=infl, expect = opt$expect) 
  out$centralBorder <- border$centralBorder
  out$outerBorder <- border$outerBorder
  out$raw.mean <-  mean(resp$x[resp$y>0])
  class(out) <- c("Para.HOF")
  out
}


"Para.HOF.list" <- function (resp, ...)   {
    out <- lapply(resp, Para, ...)
    #    class(out) <- "Para.HOF.frame"
    out
}
