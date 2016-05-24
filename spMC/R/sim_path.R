sim_path <-
function(x, data, coords, grid, radius, fixed=FALSE) {
  # Generation of conditional simulation based on path
  #
  #      x a multi_tpfit object
  #   data vector of data
  # coords coordinates matrix
  #   grid simulation points
  # radius radius to find neighbour points
  #  fixed boolean for random or fixed path algorithm

  if(!is.multi_tpfit(x)) stop("argument \"x\" must be a 'multi_tpfit' object.")

  if (missing(grid)) stop("simulation grid is missing.")
  if (missing(radius)) stop("searching radius is missing.")
  storage.mode(radius) <- "double"
  if (radius <= 0) stop("searching radius must be positive.")
  if (!is.factor(data)) data <- as.factor(data)
  if (!is.matrix(coords)) coords <- as.matrix(coords)
  if (!is.matrix(grid)) grid <- as.matrix(grid)
  nc <- dim(coords)[2]
  nr.orig <- dim(coords)[1]
  if (length(data) != nr.orig) stop("the number of data is not equal to the number of coordinates.")
  if (nc != dim(grid)[2]) stop("coordinates and simulation grid must have the same number of columns.")
  nk <- nlevels(data)
  levelLab <- levels(data)
  nrs <- dim(grid)[1] #total number of simulation
  storage.mode(coords) <- "double"
  storage.mode(grid) <- "double"

  dire.mat <- diag(, nc)
  new.coords <- coords
  new.grid <- grid
  if (!is.null(x$rotation)) {
    dire.mat <- .C('rotaxes', nc = as.integer(nc), ang = as.double(x$rotation), 
                   res = as.double(dire.mat), PACKAGE = "spMC")$res
    dire.mat <- matrix(dire.mat, nc, nc)
    new.coords <- matrix(.C('fastMatProd', nr = as.integer(nr.orig), ni = as.integer(nc),
                          mat1 = as.double(coords), nc = as.integer(nc),
                          mat2 = as.double(dire.mat), res = as.double(new.coords),
                          PACKAGE = "spMC")$res, nrow = nr.orig, ncol = nc)
    new.grid <- matrix(.C('fastMatProd', nr = as.integer(nrs), ni = as.integer(nc),
                          mat1 = as.double(grid), nc = as.integer(nc),
                          mat2 = as.double(dire.mat), res = as.double(new.grid),
                          PACKAGE = "spMC")$res, nrow = nrs, ncol = nc)
  }

  # GENERATING THE PATH TO FOLLOW
  if (fixed) {
    # Fixed Path MC Algorithm
    ord <- do.call("order", as.data.frame(grid[, nc:1]))
    path <- as.integer(as.factor(grid[, nc]))
    path <- unlist(sapply(1:max(path), function(i) if(i%%2) ord[i == path] else rev(ord[i==path])))
  }
  else {
    # Random Path MC Algorithm
    path <- sample(1:nrs)
  }
  storage.mode(path) <- "integer"

  data <- as.integer(data)
  prhat <- matrix(0, nrow = nrs, ncol = nk)
  pred <- vector("integer", nrs)

  res <- .C('pathAlg', nrs = as.integer(nrs), nrorig = as.integer(nr.orig),
            nc = as.integer(nc), coords = as.double(new.coords), gird = as.double(new.grid),
            path = as.integer(path), radius = as.double(radius), nk = as.integer(nk),
            data = as.integer(data), coefs = as.double(unlist(x$coefficients)),
            prop = as.double(x$prop), prhat = as.double(prhat), pred = as.integer(pred),
            PACKAGE = "spMC")[12:13]

  prhat <- matrix(res$prhat, nrow = nrs, ncol = nk)

  simu <- vector("integer", nrs)
  simu <- .C('tsimCate', nk = as.integer(nk), n = as.integer(nrs), prhat = as.double(prhat), 
            initSim = as.integer(simu), PACKAGE = "spMC")$initSim

  tmpfct <- 1:nk
  tmpfct <- factor(tmpfct, labels = levelLab)
  simu <- tmpfct[simu]
  pred <- tmpfct[res$pred]
  rownames(grid) <- NULL 
  rownames(pred) <- NULL
  rownames(simu) <- NULL
  rownames(prhat) <- NULL
  res <- data.frame(grid, simu, pred, prhat)
  names(res) <- c(colnames(coords), "Simulation", "Prediction", levelLab)
  tipo <- if (fixed) {"Fixed"} else {"Random"}
  attr(res, "type") <- paste(tipo, "Path Simulation")
  class(res) <- c("data.frame", "spsim")
  return(res)
}
