##' @export 
tobitw_method.lvm <- "estfun"

##' @export 
tobitw_gradient.lvm <- tobit_gradient.lvm

##' @export 
tobitw_hessian.lvm <- tobit_hessian.lvm


## tobitw2_method.lvm <- "NR"
## tobitw2_objective.lvm <- function(...) {
##   S <- tobit_gradient.lvm(...)
##   crossprod(S)[1]
## }
## tobitw2_gradient.lvm <- function(...) {
##   tobit_gradient.lvm(...)
## }
## tobitw2_hessian.lvm <- function(p,...) {
##   S <- tobit_gradient.lvm(p=p,...)
##   myfun <- function(p0) tobitw_gradient.lvm(p=p0,...)
##   H <- jacobian(myfun,p,method=lava.options()$Dmethod)
##   attributes(H)$grad <- S
##   H
## }

