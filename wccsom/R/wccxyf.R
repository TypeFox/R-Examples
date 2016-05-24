### We can later add arguments init and nhbrdist if we want
"wccxyf" <- function(data, Y, grid = somgrid(), rlen = 100,
                     alpha = c(0.05, 0.01),
                     radius = quantile(nhbrdist, 0.67),
                     xweight = 0.5, trwidth = 20, 
                     toroidal = FALSE, keep.data = TRUE)
{
  data <- as.matrix(data)
  datadims <- dimnames(data)
  dimnames(data) <- NULL
  nobj <- nrow(data)
  nvar <- ncol(data)
  nunits <- nrow(grid$pts)
  if (is.vector(Y)) Y <- matrix(Y, ncol=1)
  ny <- ncol(Y)
  predict.type <- ifelse (all(rowSums(Y) == 1),
                          "class",
                          "continuous")

  ## Fill distance matrix
  ##  if (missing(nhbrdist))
  nhbrdist <- unit.distances(grid, toroidal)
  if (toroidal) {
    if (grid$topo == "hexagonal" & (grid$ydim %% 2 == 1))
      stop("Error: uneven number of rows (y) in hexagonal toroidal grid")
    radius <- radius*0.5
  }

  ## Set weights for crosscorrelation function
  if (trwidth > 0)
    wghts <- 1 - (0:trwidth)/trwidth
  else
    wghts <- 1

  ## Initialise codebook vectors
  ## In cases where smaller classes are swamped by larger ones we
  ## maybe should augment the small classes a little bit... TODO
  starters <- sample(nobj, nunits, replace=FALSE)
  codes <- data[starters,]
  codeYs <- Y[starters,]
  ## In the C-code, it is easier to have data
  ## and codes transposed.
  data <- t(data)
  codes <- t(codes)
  codeYs <- t(codeYs)
  Y <- t(Y)
  
  ## declare space for arrays
  changes <- rep(0, rlen*2)
  xdists <- ydists <- rep(0, nunits)  
  acors <- rep(0, nunits)
  data.acors <- rep(0, nobj)

  if (predict.type == "continuous") {
    res <- .C("WCCXYF_Eucl",
              data = as.double(data),
              Ys = as.double(Y),
              codes = as.double(codes),
              codeYs = as.double(codeYs),
              nhbrdist = as.double(nhbrdist),
              alpha = as.double(alpha),
              radius = as.double(radius),
              xweight = as.double(xweight),
              trwdth = as.integer(trwidth),
              wghts = as.double(wghts),
              data.acors = as.double(data.acors),
              acors = as.double(acors),
              changes = as.double(changes),
              xdists = as.double(xdists),
              ydists = as.double(ydists),
              n = as.integer(nobj),
              px = as.integer(nvar),
              py = as.integer(ny),
              ncodes = as.integer(nunits),
              rlen = as.integer(rlen),
              PACKAGE = "wccsom")
  } else {
    res <- .C("WCCXYF_Tani",
              data = as.double(data),
              Ys = as.double(Y),
              codes = as.double(codes),
              codeYs = as.double(codeYs),
              nhbrdist = as.double(nhbrdist),
              alpha = as.double(alpha),
              radius = as.double(radius),
              xweight = as.double(xweight),
              trwdth = as.integer(trwidth),
              wghts = as.double(wghts),
              data.acors = as.double(data.acors),
              acors = as.double(acors),
              changes = as.double(changes),
              xdists = as.double(xdists),
              ydists = as.double(ydists),
              n = as.integer(nobj),
              px = as.integer(nvar),
              py = as.integer(ny),
              ncodes = as.integer(nunits),
              rlen = as.integer(rlen),
              PACKAGE = "wccsom")
  }

  changes <- matrix(res$changes, ncol=2)
  codes <- matrix(res$codes, nvar, nunits)
  codeYs <- matrix(res$codeYs, ny, nunits)
  rownames(codeYs) <- rownames(Y)
  acors <- res$acors
  data.acors <- res$data.acors

  if (keep.data) {
    classif <- wccs <- rep(0, nobj)
    res <-  .C("wccassign",
               data = as.double(data),
               data.acors = as.double(data.acors),
               codes = as.double(codes),
               acors = as.double(acors),
               as.integer(trwidth),
               wghts = as.double(wghts),
               classif = as.integer(classif),
               wccs = as.double(wccs),
               as.integer(nobj),
               as.integer(nvar),
               as.integer(nunits),
               PACKAGE = "wccsom")
    unit.classif <- res$classif
    wccs <- res$wccs
    
    structure(list(data = t(data), Y = t(Y), predict.type = predict.type,
                   grid = grid, changes = changes,
                   codes = t(codes), codeYs = t(codeYs),
                   trwdth = trwidth,
                   unit.classif = unit.classif, wccs = wccs,
                   data.acors = data.acors, acors = acors,
                   toroidal = toroidal, xweight = xweight),
              class = c("wccsom", "SOM"))
  } else {
    structure(list(predict.type = predict.type,
                   grid = grid, changes = changes,
                   codes = t(codes), codeYs = t(codeYs),
                   trwdth = trwidth, acors = acors,
                   toroidal = toroidal, xweight = xweight), 
              class = c("wccsom", "SOM"))
  }
}

