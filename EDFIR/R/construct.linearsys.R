construct.linearsys <-
function(predatorframe, preyframe) {
  preyhull = build.normalhull(preyframe)
  if (dim(predatorframe)[1] > 1) {
    predatorhull = build.hull.list(predatorframe)
  } else {
    predator = predatorframe[1,]
  }
  ## Predator values on the convex hull
  n = dim(predator)[1]
  ## Now create the linear system for the input to the program
  nrows = dim(preyhull)[1]
  ncols = dim(preyhull)[2]
  ## Create the "A" of Ax <= b. The "A" in this case needs to
  ## be doubled up since x < 0, and repeated for every predator
  ## ocurring in the hull.
  A = preyhull[,1:(ncols-1)]
  inA = A[rep(1:nrows, n),]
  b = -preyhull[,ncols]
  inb = apply(predator, 1, function(x) { b - A%*%x })
  ## we take the negative since the normals are outward facing (see the Qhull documentation
  ## for output with the "n" flag for more information)
  list(A = inA, b = as.vector(inb))
}
