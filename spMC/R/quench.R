quench <-
function(x, data, coords, sim, GA = FALSE, optype = c("param", "fullprobs", "semiprobs", "coordprobs"), max.it = 1000, knn = 12) {
  # Perform quenching algorithm for reducing artifacts on simulations
  #
  #        x a multi_tpfit object
  #     data vector of data
  #   coords coordinates matrix
  #      sim simulation results (tipically an "spsim" object)
  #       GA boolean (if TRUE genetic algorithm is applied rather than simulated annealing)
  #   optype character with the objective function to minimize after the simulation
  #   max.it maximum number of iteration for the optimization method
  #      knn number of k-nearest neighbours

  if(!is.multi_tpfit(x)) stop("argument \"x\" must be a 'multi_tpfit' object.")
  if(!all(class(sim) %in% c("data.frame", "spsim"))) stop("argument \"sim\" must be an 'spsim' object.")
  if (!is.matrix(coords)) coords <- as.matrix(coords)
  nc <- dim(coords)[2]
  nr <- dim(coords)[1]
  grid <- sim[, 1L:nc]
  if (!is.factor(data)) data <- as.factor(data)
  if (!is.matrix(grid)) grid <- as.matrix(grid)

  probs <- sim[, -1L:-(nc + 2L)]

  nrs <- dim(grid)[1]
  nk <- nlevels(data)
  if (is.null(knn)) {knn <- nr} else {knn <- as.integer(knn)}
  levLabels <- levels(data)
  data <- as.integer(data)
  storage.mode(coords) <- "double"
  storage.mode(grid) <- "double"
  if (length(optype) > 1) optype <- optype[1]

  initSim <- sim$Simulation
  pred <- sim$Prediction

  dire.mat <- diag(, nc)
  new.coords <- coords
  new.grid <- grid
  if (!is.null(x$rotation)) {
    dire.mat <- .C('rotaxes', nc = as.integer(nc), ang = as.double(x$rotation), 
                   res = as.double(dire.mat), PACKAGE = "spMC")$res
    dire.mat <- matrix(dire.mat, nc, nc)
    new.coords <- matrix(.C('fastMatProd', nr = as.integer(nr), ni = as.integer(nc),
                          mat1 = as.double(coords), nc = as.integer(nc),
                          mat2 = as.double(dire.mat), res = as.double(new.coords),
                          PACKAGE = "spMC")$res, nrow = nr, ncol = nc)
    new.grid <- matrix(.C('fastMatProd', nr = as.integer(nrs), ni = as.integer(nc),
                          mat1 = as.double(grid), nc = as.integer(nc),
                          mat2 = as.double(dire.mat), res = as.double(new.grid),
                          PACKAGE = "spMC")$res, nrow = nrs, ncol = nc)
  }

  # FINDING THE k-NEAREST NEIGHBOURS #
  indices <- matrix(0L, nrow = knn, ncol = nrs)
  indices <- matrix(.C('knear', nc = as.integer(nc), nr = as.integer(nr),
                       coords = as.double(new.coords), nrs = as.integer(nrs),
                       grid = as.double(new.grid), knn = as.integer(knn),
                       indices = as.integer(indices), PACKAGE = "spMC")$indices,
                    nrow = knn, ncol = nrs)

  # SORTING SIMULATION GRID #
  path <- do.call("order", as.data.frame(t(indices)))
  indices <- indices[, path]
  grid <- grid[path, ]
  new.grid <- new.grid[path, ]
  groups <- !duplicated(t(indices))
  groups <- cumsum(groups)
  initSim <- initSim[path]

  # OPTIMIZATION PROCEDURE #
  toOptim <- function(x, mySim, grid) {
    xnew <- list()
    xnew$coefficients <- lapply(1:nc, function(i) {
      ml <- mlen(mySim, grid, loc.id[, i], dire.mat[i, ])
      Rmat <- embed_MC(mySim, grid, loc.id[, i], dire.mat[i, ])
      diag(Rmat) <- -1
      Rmat <- diag(1 / ml) %*% Rmat
      return(Rmat)
    })
    xnew$prop <- table(mySim)
    xnew$prop <- xnew$prop / sum(xnew$prop)
    if(length(unlist(x$coefficients)) != length(unlist(xnew$coefficients))) return(Inf)
    if(length(unlist(x$prop)) != length(unlist(xnew$prop))) return(Inf)
    return(sum((unlist(xnew$coefficients) - unlist(x$coefficients))^2) + sum((xnew$prop - x$prop)^2))
  }
  if(optype == "fullprobs") {
    toOptim <- function(x, mySim, grid) {
      res <- 0
      res <- .C('objfun', nrs = as.integer(nrs), nk = as.integer(nk), nc = as.integer(nc),
         mySim = as.integer(mySim), grid = as.double(grid),
         coef = as.double(unlist(x$coefficients)), prop = as.double(x$prop),
         res = as.double(res), PACKAGE = "spMC")$res
      return(res)
    }
  }
  if(optype == "semiprobs") {
    sknn <- knn + 1L
    storage.mode(sknn) <- "integer"
    if(nrs < sknn) sknn <- nrs
    indicesim <- matrix(0L, nrow = sknn, ncol = nrs)
    indicesim <- matrix(.C('knear', nc = as.integer(nc), nr = as.integer(nrs),
                         coords = as.double(new.grid), nrs = as.integer(nrs),
                         grid = as.double(new.grid), knn = as.integer(sknn),
                         indices = as.integer(indicesim),
                         PACKAGE = "spMC")$indices,
                      nrow = sknn, ncol = nrs)
    sknn <- sknn - 1L
    toOptim <- function(x, mySim, grid) {
      res <- 0
      res <- .C('fastobjfun', knn = as.integer(sknn),
         indices = as.integer(indicesim[-1L, ]), nrs = as.integer(nrs), nk = as.integer(nk),
         nc = as.integer(nc), nr = as.integer(nrs), mySim = as.integer(mySim),
         grid = as.double(grid), coef = as.double(unlist(x$coefficients)),
         prop = as.double(x$prop), data = as.integer(mySim), coords = as.double(grid),
         res = as.double(res), PACKAGE = "spMC")$res
      return(res)
    }
  }
  if(optype == "coordprobs") {
    toOptim <- function(x, mySim, grid) {
      res <- 0
      res <- .C('fastobjfun', knn = as.integer(knn), indices = as.integer(indices),
         nrs = as.integer(nrs), nk = as.integer(nk), nc = as.integer(nc), nr = as.integer(nr),
         mySim = as.integer(mySim), grid = as.double(grid),
         coef = as.double(unlist(x$coefficients)), prop = as.double(x$prop),
         data = as.integer(data), coords = as.double(coords), res = as.double(res),
         PACKAGE = "spMC")$res
      return(res)
    }
  }
  if (max.it > 0 & !(optype %in% c("fullprobs", "semiprobs", "coordprobs"))) {
    loc.id <- apply(dire.mat, 1, function(d) which_lines(grid, d, tolerance = x$tolerance))
  }
  storage.mode(initSim) <- "integer"
  storage.mode(max.it) <- "integer"
  Rnv <- new.env()
  Rnv <- parent.env(Rnv)

  if (!GA) { # SIMULATED ANNEALING #
    old <- .Call("annealingSIM", max.it, initSim, x, grid, quote(toOptim(x, pp, grid)), Rnv, PACKAGE = "spMC")
  }
  else {     # GENETIC ALGORITHM #
    old <- .Call("geneticSIM", max.it, initSim, x, grid, quote(toOptim(x, pp, grid)), Rnv, PACKAGE = "spMC")
  }
  tmpfct <- 1:nk
  tmpfct <- factor(tmpfct, labels = levLabels)
  old <- tmpfct[old]
  res <- data.frame(grid, old, pred, probs)
  names(res) <- c(colnames(coords), "Simulation", "Prediction", levLabels)
  res[path, ] <- res
  attr(res, "type") <- attr(sim, "type")
  class(res) <- c("data.frame", "spsim")
  return(res)
}
