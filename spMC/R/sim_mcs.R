sim_mcs <-
function(x, data, coords, grid, knn = NULL) {
  # Generation of Multinomial Categorical Simulation (MCS)
  #
  #      x a multi_tpfit object
  #   data vector of data
  # coords coordinates matrix
  #   grid simulation points
  #    knn number of k-nearest neighbours (if NULL all data are neighbours)
  # radius radius to find neighbour points

  if(!is.multi_tpfit(x)) stop("argument \"x\" must be a 'multi_tpfit' object.")

  if (missing(grid)) stop("simulation grid is missing.")
  if (!is.factor(data)) data <- as.factor(data)
  if (!is.matrix(coords)) coords <- as.matrix(coords)
  if (!is.matrix(grid)) grid <- as.matrix(grid)
  storage.mode(coords) <- "double"
  storage.mode(grid) <- "double"

  hmany <- length(data)
  nc <- dim(coords)[2]
  nr.orig <- dim(coords)[1]
  if (hmany != nr.orig) stop("the number of data is not equal to the number of coordinates")
  if (nc != dim(grid)[2]) stop("coordinates and simulation grid must have the same number of columns")
  if (!is.null(knn)) {
    knn <- if ((knn > 0) && (knn < hmany)) { as.integer(ceiling(knn)) } else { NULL }
  }
  nk <- nlevels(data)
  levelLab <- levels(data)
  nrs <- dim(grid)[1] #total number of simulation

  data <- as.integer(data)
  prhat <- matrix(x$prop, nk, nrs)
  dire.mat <- diag(, nc)
  if (!is.null(x$rotation)) {
    dire.mat <- .C('rotaxes', nc = as.integer(nc), ang = as.double(x$rotation),
                   res = as.double(dire.mat), PACKAGE = "spMC")$res
    dire.mat <- matrix(dire.mat, nc, nc)
  }
  if (is.null(knn)) {
    path <- 1:nrs
    # APPROXIMATING PROBABILITIES WITH ALL DATA #
    prhat <- .C('jointProbsMCS', coords = as.double(coords), hmany = as.integer(hmany),
                grid = as.double(grid), nrs = as.integer(nrs), nc = as.integer(nc),
                nk = as.integer(nk), ndata = as.integer(data),
                coefs = as.double(unlist(x$coefficients)), matdir = as.double(dire.mat),
                rota = as.integer(!is.null(x$rotation)), pProbs = as.double(prhat), 
                PACKAGE = "spMC")$pProbs
  }
  else {
    new.coords <- coords
    new.grid <- grid
    if (!is.null(x$rotation)) {
      new.coords <- matrix(.C('fastMatProd', nr = as.integer(nr.orig), ni = as.integer(nc),
                              mat1 = as.double(coords), nc = as.integer(nc),
                              mat2 = as.double(dire.mat), res = as.double(new.coords),
                              PACKAGE = "spMC")$res, nrow = nr.orig, ncol = nc)
      new.grid <- matrix(.C('fastMatProd', nr = as.integer(nrs), ni = as.integer(nc),
                            mat1 = as.double(grid), nc = as.integer(nc),
                            mat2 = as.double(dire.mat), res = as.double(new.grid),
                            PACKAGE = "spMC")$res, nrow = nrs, ncol = nc)
    }
    # FINDING THE k-NEAREST NEIGHBOURS #
    indices <- matrix(0L, nrow = knn, ncol = nrs)
    indices <- matrix(.C('knear', nc = as.integer(nc), nr = as.integer(nr.orig),
                         coords = as.double(new.coords), nrs = as.integer(nrs),
                         grid = as.double(new.grid), knn = as.integer(knn),
                         indices = as.integer(indices), PACKAGE = "spMC")$indices,
                      nrow = knn, ncol = nrs)
    # SORTING SIMULATION GRID #
    path <- do.call("order", as.data.frame(t(indices)))
    indices <- indices[, path]
    grid <- grid[path, ]
    # APPROXIMATING PROBABILITIES WITH k-NEAREST NEIGHBOURS #
    prhat <- .C('KjointProbsMCS', coords = as.double(new.coords), hmany = as.integer(hmany),
                grid = as.double(new.grid), nrs = as.integer(nrs), nc = as.integer(nc),
                nk = as.integer(nk), ndata = as.integer(data), knn = as.integer(knn),
                coefs = as.double(unlist(x$coefficients)), indices = as.integer(indices),
                pProbs = as.double(prhat), PACKAGE = "spMC")$pProbs
  }
  prhat <- matrix(prhat, nrs, nk, byrow = TRUE)
  sim <- vector("integer", nrs)

  sim <- .C('tsimCate', nk = as.integer(nk), n = as.integer(nrs), prhat = as.double(prhat), 
            initSim = as.integer(sim), PACKAGE = "spMC")$initSim
  pred <- apply(prhat, 1, which.max)

  rownames(grid) <- NULL
  rownames(sim) <- NULL
  rownames(pred) <- NULL
  rownames(prhat) <- NULL
  tmpfct <- 1:nk
  tmpfct <- factor(tmpfct, labels = levelLab)
  sim <- tmpfct[sim]
  pred <- tmpfct[pred]
  res <- data.frame(grid, sim, pred, prhat)
  names(res) <- c(colnames(coords), "Simulation", "Prediction", levelLab)
  res[path, ] <- res
  attr(res, "type") <- "Multinomial Categorical Simulation"
  class(res) <- c("data.frame", "spsim")
  return(res)
}
