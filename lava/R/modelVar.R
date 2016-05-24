
###{{{ modelVar

##' @export
`modelVar` <-
  function(x,p,...) UseMethod("modelVar")

##' @export
modelVar.lvmfit <- function(x, p=pars(x), ...) modelVar(Model(x),p=p,...)

##' @export
modelVar.lvm <- function(x,p,data,...) {
  pp <- modelPar(x,p)
  res <- moments(x, p=p, data=data,...)
  attr(res, "pars") <- pp$p
  attr(res, "meanpar") <- pp$meanpar
  attr(res, "epar") <- pp$epar
  res
}
###}}} modelVar
