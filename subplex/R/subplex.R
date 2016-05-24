subplex <- function (par, fn, control = list(), hessian = FALSE, ...) {
## set default control parameters
  con <- list(
              reltol = .Machine$double.eps,
              maxit = 10000,
              parscale = rep.int(1,length(par))
              )
  namc <- names(control)[names(control)%in%names(con)]
  con[namc] <- control[namc]
  .Call(
        call_subplex,
        par,
        match.fun(fn),
        tol=con$reltol,
        maxnfe=con$maxit,
        scale=con$parscale,
        hessian,
        environment(fn),
        pairlist(...)
        )
}
