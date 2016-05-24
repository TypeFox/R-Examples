optimiDefaultOptions <-
function() {
  ## options: trace, maximum iteration, fnscale, reltol, default optimiser
  ## return (list(trace=TRUE, maxit=1000, fnscale=1e1, reltol=1e-4, optimiser="CG", gradcheck=FALSE, hessian=FALSE))
  return (list(maxit=3000, ln=c(0,2), xtol=1e-4, fnTol=1e-4, optimiser="SCG", gradcheck=FALSE, display=TRUE))
}
