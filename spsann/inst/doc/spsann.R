## ----eval=FALSE----------------------------------------------------------
#  # OPTIMIZING A SAMPLE CONFIGURATION IN FOUR STEPS
#  # Step 1. Load and pre-process the data
#  # ...
#  # Step 2. Set control parameters
#  # ...
#  # Step 3. Execute the simulated annealing algorithm
#  # ...
#  # ... Be prepared to wait!!!
#  # ...
#  # Step 4. Evaluate the optimized sample configuration
#  # ...

## ----eval=FALSE----------------------------------------------------------
#  # Load and pre-process the data
#  data(meuse.grid, package = "sp")
#  boundary <- meuse.grid
#  sp::coordinates(boundary) <- c("x", "y")
#  sp::gridded(boundary) <- TRUE
#  boundary <- rgeos::gUnaryUnion(as(boundary, "SpatialPolygons"))
#  candi <- meuse.grid[, 1:2]
#  covars <- meuse.grid[, 6:7]
#  
#  # Set control parameters
#  schedule <- scheduleSPSANN(initial.temperature = 0.5)
#  set.seed(2001)
#  
#  # Execute the simulated annealing algorithm
#  res <- optimDIST(
#    points = 30, candi = candi, covars = covars, use.coords = TRUE,
#    schedule = schedule, plotit = TRUE, boundary = boundary)
#  
#  # Evaluate the optimized sample configuration
#  objSPSANN(res)
#  objDIST(
#    points = res, candi = candi, covars = covars, use.coords = TRUE)
#  plot(res, boundary = boundary)

## ----eval=FALSE----------------------------------------------------------
#  # Load and pre-process the data
#  data(meuse.grid, package = "sp")
#  boundary <- meuse.grid
#  sp::coordinates(boundary) <- c("x", "y")
#  sp::gridded(boundary) <- TRUE
#  boundary <- rgeos::gUnaryUnion(as(boundary, "SpatialPolygons"))
#  candi <- meuse.grid[, 1:2]
#  
#  # Set control parameters
#  schedule <- scheduleSPSANN(initial.temperature = 500)
#  set.seed(2001)
#  
#  # Execute the simulated annealing algorithm
#  res <- optimPPL(
#    points = 30, candi = candi, pairs = TRUE, schedule = schedule,
#    plotit = TRUE, boundary = boundary)
#  
#  # Evaluate the optimized sample configuration
#  objSPSANN(res)
#  objPPL(points = res, pairs = TRUE, candi = candi)
#  countPPL(points = res, candi = candi, pairs = TRUE)
#  plot(res, boundary = boundary)

## ----eval=FALSE----------------------------------------------------------
#  # Load and pre-process the data
#  data(meuse.grid, package = "sp")
#  boundary <- meuse.grid
#  sp::coordinates(boundary) <- c("x", "y")
#  sp::gridded(boundary) <- TRUE
#  boundary <- rgeos::gUnaryUnion(as(boundary, "SpatialPolygons"))
#  candi <- meuse.grid[, 1:2]
#  
#  # Set control parameters
#  schedule <- scheduleSPSANN(
#    initial.acceptance = 0, initial.temperature = 0.01)
#  set.seed(2001)
#  
#  # Execute the simulated annealing algorithm
#  res <- optimMSSD(
#    points = 30, candi = candi, schedule = schedule, plotit = TRUE,
#    boundary = boundary)
#  
#  # Evaluate the optimized sample configuration
#  objSPSANN(res)
#  objMSSD(candi = candi, points = res)
#  plot(res, boundary = boundary)

## ----eval=FALSE----------------------------------------------------------
#  # Load and pre-process the data
#  data(meuse.grid, package = "sp")
#  boundary <- meuse.grid
#  sp::coordinates(boundary) <- c("x", "y")
#  sp::gridded(boundary) <- TRUE
#  boundary <- rgeos::gUnaryUnion(as(boundary, "SpatialPolygons"))
#  candi <- meuse.grid[, 1:2]
#  covars <- as.data.frame(meuse.grid)
#  
#  # Set control parameters
#  vgm <- gstat::vgm(
#    psill = 10, model = "Exp", range = 500, nugget = 8)
#  schedule <- scheduleSPSANN(initial.temperature = 10)
#  set.seed(2001)
#  
#  # Execute the simulated annealing algorithm
#  res <- optimMKV(
#    points = 30, candi = candi, covars = covars, vgm = vgm,
#    eqn = z ~ soil, plotit = TRUE, boundary = boundary,
#    schedule = schedule)
#  
#  # Evaluate the optimized sample configuration
#  objSPSANN(res)
#  objMKV(
#    points = res, candi = candi, covars = covars,
#    eqn = z ~ soil, vgm = vgm)
#  plot(res, boundary = boundary)

## ----eval=FALSE----------------------------------------------------------
#  # Load and pre-process the data
#  data(meuse.grid, package = "sp")
#  boundary <- meuse.grid
#  sp::coordinates(boundary) <- c("x", "y")
#  sp::gridded(boundary) <- TRUE
#  boundary <- rgeos::gUnaryUnion(as(boundary, "SpatialPolygons"))
#  candi <- meuse.grid[, 1:2]
#  
#  # Set control parameters
#  schedule <- scheduleSPSANN(
#    initial.temperature = 30, x.max = 1540, y.max = 2060,
#    x.min = 0, y.min = 0, cellsize = 40)
#  objUSER <- function (points, lags, n_lags, n_pts) {
#    dm <- SpatialTools::dist1(points[, 2:3])
#    ppl <- vector()
#    for (i in 1:n_lags) {
#      n <- which(dm > lags[i] & dm <= lags[i + 1], arr.ind = TRUE)
#      ppl[i] <- length(unique(c(n)))
#    }
#    distri <- rep(n_pts, n_lags)
#    res <- sum(distri - ppl)
#  }
#  lags <- seq(1, 1000, length.out = 10)
#  set.seed(2001)
#  
#  # Execute the simulated annealing algorithm
#  res <- optimUSER(
#    points = 30, fun = objUSER, lags = lags, n_lags = 9,
#    n_pts = 10, candi = candi, schedule = schedule,
#    plotit = TRUE, boundary = boundary)
#  
#  # Evaluate the optimized sample configuration
#  objSPSANN(res)
#  countPPL(res, candi = candi, lags = lags)
#  plot(res, boundary = boundary)

