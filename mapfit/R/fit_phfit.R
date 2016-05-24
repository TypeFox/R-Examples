

phfit.point <- function(ph, x, weights, method = c("all", "increment"),
  lbound = 1, ubound = NULL, control = list(), verbose = list(), ...) {
  data <- phfit.time.data.frame(time=x, weights=weights)
  switch(class(ph),
    "ph"=phfit.gen(ph=ph, data=data, control=control, verbose=verbose, ...),
    "cf1"=phfit.cf1(ph=ph, data=data, control=control, verbose=verbose, ...),
    "herlang"={
      phsize <- sum(ph@shape)
      if (is.null(ubound)) {
        ubound <- phsize
      }
      phfit.herlang(phsize=phsize, data=data, method=method, lbound=lbound, ubound=ubound,
        control=control, verbose=verbose, ...)
    })
}

phfit.group <- function(ph, counts, breaks, intervals, instant,
 method = c("all", "increment"), lbound = 1, ubound = NULL, control = list(), verbose = list(), ...) {
  data <- phfit.group.data.frame(counts=counts, breaks=breaks, difftime=intervals, instant=instant)
  switch(class(ph),
    "ph"=phfit.gen(ph=ph, data=data, control=control, verbose=verbose, ...),
    "cf1"=phfit.cf1(ph=ph, data=data, control=control, verbose=verbose, ...),
    "herlang"={
      phsize <- sum(ph@shape)
      if (is.null(ubound)) {
        ubound <- phsize
      }
    	phfit.herlang(phsize=phsize, data=data, method=method, lbound=lbound, ubound=ubound,
    		control=control, verbose=verbose, ...)
    })
}

phfit.density <- function(ph, f, deformula = zero.to.inf, weight.zero = .Machine$double.eps,
  weight.reltol = sqrt(.Machine$double.eps),
  method = c("all", "increment"), lbound = 1, ubound = NULL, 
  control = list(), verbose = list(), ...) {
  x <- deformula.weight(f, deformula, weight.zero, weight.reltol, ...)
  data <- phfit.time.data.frame(time=x$x, weights=x$weights)
  switch(class(ph),
    "ph"=phfit.gen(ph=ph, data=data, control=control, verbose=verbose, ...),
    "cf1"=phfit.cf1(ph=ph, data=data, control=control, verbose=verbose, ...),
    "herlang"={
      phsize <- sum(ph@shape)
      if (is.null(ubound)) {
        ubound <- phsize
      }
      phfit.herlang(phsize=phsize, data=data, method=method, lbound=lbound, ubound=ubound,
        control=control, verbose=verbose, ...)
    })
}
