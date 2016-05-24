sab.init <- function(constr,
    thin.fn = function(n) { ceiling(log(n + 1)/4 * n^3) },
    thin = NULL,
    x0.randomize = FALSE, x0.method="slacklp",
    x0 = NULL,
    eliminate=TRUE) {
  state <- har.init(constr, thin.fn, thin, x0.randomize, x0.method, x0, eliminate)
  state$i0 <- findFace(state$x0, state$constr) # find the closest face of the polytope
  state
}

sab.run <- function(state, n.samples) {
  result <- with(state, {
    n <- length(x0) - 1
    if (n == 0) {
      list(samples = matrix(rep(basis$translate, each=n.samples), nrow=n.samples), xN = 1, iN = 1)
    } else {
      sab(x0, i0, constr, N=n.samples * thin, thin=thin, homogeneous=TRUE, transform=transform)
    }
  })
  state$x0 <- result$xN
  state$i0 <- result$iN
  list(state = state, samples = result$samples, faces = result$faces)
}

shakeandbake <- function(constr,
    n.samples=1E4,
    thin.fn = function(n) { ceiling(log(n + 1)/4 * n^3) },
    thin = NULL,
    x0.randomize = FALSE, x0.method="slacklp",
    x0 = NULL,
    eliminate=TRUE) {
  state <- sab.init(constr, thin.fn, thin, x0.randomize, x0.method, x0, eliminate)
  result <- sab.run(state, n.samples)
  result$samples
}
