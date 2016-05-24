"sim_ik" <- 
function(x, data, coords, grid, knn = 12, ordinary = TRUE) {
  # Generation of conditional simulation based on Kriging
  #
  #        x a multi_tpfit object
  #     data vector of data
  #   coords coordinates matrix
  #     grid simulation points
  #      knn number of k-nearest neighbours
  # ordinary boolean (if TRUE ordinary Kriging is applied rather than simple Kriging)
  #       GA boolean (if TRUE genetic algorithm is applied rather than simulated annealing)
  #   optype character with the objective function to minimize after the simulation
  #   max.it maximum number of iteration for the optimization method

  ordinary <- as.logical(ordinary)
  if(!is.multi_tpfit(x)) stop("argument \"x\" must be a 'multi_tpfit' object.")

  if (missing(grid)) stop("simulation grid is missing.")
  if (!is.factor(data)) data <- as.factor(data)
  if (!is.matrix(coords)) coords <- as.matrix(coords)
  if (!is.matrix(grid)) grid <- as.matrix(grid)
  nc <- dim(coords)[2]
  nr <- dim(coords)[1]
  if (is.null(knn)) {knn <- nr} else {knn <- as.integer(knn)}
  if (length(data) != nr) stop("the number of data is not equal to the number of coordinates")
  if (nc != dim(grid)[2]) stop("coordinates and simulation grid must have the same number of columns")
  if (nr < knn) knn <- nr
  nrs <- dim(grid)[1]
  nk <- nlevels(data)
  levLabels <- levels(data)
  data <- as.integer(data)
  storage.mode(coords) <- "double"
  storage.mode(grid) <- "double"
#   if (length(optype) > 1) optype <- optype[1]

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

  # KRIGING PROCEDURE #
  probs <- matrix(0, nrow = nrs, ncol = nk)
  probs <- matrix(.C('getIKPrbs', ordinary = as.integer(ordinary),
                     indices = as.integer(indices), groups = as.integer(groups), 
                     knn = as.integer(knn), nc = as.integer(nc), nr = as.integer(nr),
                     nrs = as.integer(nrs), data = as.integer(data),
                     coords = as.double(new.coords), grid = as.double(new.grid),
                     nk = as.integer(nk), coef = as.double(unlist(x$coefficients)),
                     prop = as.double(x$prop), probs = as.double(probs),
                     PACKAGE = "spMC")$probs,
                  nrow = nrs, ncol = nk)

  # PREDICTION AND SIMULATION PROCEDURE #
  pred <- apply(probs, 1, which.max)
  initSim <- vector("integer", nrs)
  initSim <- .C('tsimCate', nk = as.integer(nk), n = as.integer(nrs), prhat = as.double(probs),
     initSim = as.integer(initSim), PACKAGE = "spMC")$initSim

  tmpfct <- 1:nk
  tmpfct <- factor(tmpfct, labels = levLabels)
  old <- tmpfct[initSim]
  pred <- tmpfct[pred]
  res <- data.frame(grid, old, pred, probs)
  names(res) <- c(colnames(coords), "Simulation", "Prediction", levLabels)
  res[path, ] <- res
  attr(res, "type") <- paste(ifelse(ordinary, "Ordinary", "Simple"), 
                             "Indicator Kriging Simulation")
  class(res) <- c("data.frame", "spsim")
  return(res)
}
