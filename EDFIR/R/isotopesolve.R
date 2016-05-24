isotopesolve <-
function(predatorframe, preyframe, maximize = "min") {
  out = construct.linearsys(predatorframe, preyframe)
  d = dim(predatorframe)[2]
  ## x has the form
  ## x = [x+, x-], so the value of a single coordinate of the shift
  ## vector will be something like (x[i] - x[i + dim]), where
  ## dim is the number of isotopes (i.e. dimensions).
  ## for the negative variables
  A = out$A
  A1 = cbind(A, -A)
  lpobj = lp(direction = maximize, objective.in = rep(1,2*d), const.mat = A1, const.dir = rep("<=", length(out$b)),
    const.rhs = out$b)
  lpobj
}
