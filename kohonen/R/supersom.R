"supersom" <- function(data,
                       grid = somgrid(),
                       rlen = 100,
                       alpha = c(0.05, 0.01),
                       radius = quantile(nhbrdist, 0.67) * c(1, -1),
                       contin,
                       toroidal = FALSE,
                       n.hood,
                       whatmap = NULL,
                       weights = 1,
                       maxNA.fraction = .5,
                       keep.data = TRUE)
{
  if (length(weights) == 1) weights <- rep(weights, length(data))

  whatmap <- check.whatmap(data, whatmap)
  whatmap <- whatmap[weights[whatmap] != 0]
  nmat <- length(whatmap)

  orig.data <- data
  data <- data[whatmap]
  weights <- weights[whatmap]
  orig.weights <- rep(0, length(orig.data))
  orig.weights[whatmap] <- weights

  if (length(weights) != length(data))
    stop("number of weights should equal the number of data matrices")

  if (abs(sum(weights)) < .Machine$double.eps) {
    stop("weights sum to zero")
  } else {
    weights <- weights / sum(weights)
  }

  if (!is.list(data) | !all(sapply(data, is.matrix) | sapply(data, is.factor)))
    stop("data should be a list of data matrices or factors")
  if (any(sapply(data, is.factor))) {
    data[sapply(data, is.factor)] <-
      lapply(data[sapply(data, is.factor)], classvec2classmat)
  }

  if (!all(sapply(data, is.numeric)))
    stop("Argument data should be numeric")

  if (missing(contin)) {
    ## contin == FALSE if this is a class label or membership
    ## estimate; in that case variables can be interpreted as fractions,
    ## i.e. they should not contain NAs and sum to 1
    contin <- rep(NA, length(orig.data))
    names(contin) <- names(orig.data)
    contin[whatmap] <- sapply(data,
                              function(x) {
                                !any(is.na(x)) &&
                                  any(abs(rowSums(x) - 1) > 1e-8)})
  } else {
    if (length(contin) == 1)
      contin <- rep(contin, length(orig.data))
    
    if (length(contin) != length(orig.data))
      stop("incorrect length of contin parameter")
  }
    
  ## remove NAs: individual NAs are allowed but rows or columns
  ## containing more than maxNA.fraction of NAs are removed. Columns
  ## are only removed in one data matrix.
  nacols <- lapply(data,
                   function(x)
                   which(apply(x, 2, function(y)
                               (sum(is.na(y)) / length(y)) > maxNA.fraction)))
  for (i in 1:nmat)
    if (length(nacols[[i]]) != 0) {
      warning(paste("removing", length(nacols[[i]]),
                    "NA columns from training data matrix", i, "\n"))
      data[[i]] <- data[[i]][,-nacols[[i]], drop=FALSE]
    }
  
  ## rows of NAs result in the object being removed from the data set.
  narows <- lapply(data,
                   function(x)
                   which(apply(x, 1,
                               function(y)
                               (sum(is.na(y)) / length(y)) > maxNA.fraction)))
  narows <- unique(unlist(narows))
  if (length(narows) > 0) {
    warning(paste("removing", length(narows),
                  "NA objects from the training data\n"))
    for (i in 1:nmat)
      data[[i]] <- data[[i]][-narows, , drop=FALSE]
  } else {
    narows <- NULL
  }

  ## for (i in 1:nmat)
  ##   dimnames(data[[i]]) <- NULL
  nobjects <- unique(sapply(data, nrow))
  if (length(nobjects) > 1)
    stop("unequal numbers of objects in data list")
  nvar <- sapply(data, ncol)

  nNA <- sapply(data, function(x) apply(x, 1, function(y) sum(is.na(y))))
  
  if (missing(n.hood)) {
    n.hood <- switch(grid$topo,
                     hexagonal = "circular",
                     rectangular = "square")
  } else {
    n.hood <- match.arg(n.hood, c("circular", "square"))
  }
  grid$n.hood <- n.hood
  ng <- nrow(grid$pts)
  nhbrdist <- unit.distances(grid, toroidal)
  ## temporary variables for storing distances of an object to all units
  unitdistances <- matrix(0, ng, nmat)
  if (length(radius) == 1) radius <- sort(radius * c(1, -1), decreasing = TRUE)
  
  ## initialisation
  starters <- sample(1:nobjects, ng, replace = FALSE)
  init <- vector("list", nmat)
  for (i in 1:nmat) {
    init[[i]] <- data[[i]][starters,]
    nastarters <- which(is.na(init[[i]])) # fill in NAs with random numbers
    init[[i]][nastarters] <- rnorm(length(nastarters))
  }
  
  changes <- rep(0, rlen*nmat)
  maxdists <- rep(0, nmat)
  
  ## change data and init into large concatenated matrices
  data2 <- matrix(unlist(data), nrow=nobjects)
  init <- matrix(unlist(init), nrow=ng)
  
  res <- .C("supersom",
            data = as.double(data2),
            codes = as.double(init),
            nhbrdist = as.double(nhbrdist),
            alpha = as.double(alpha),
            radii = as.double(radius),
            weights = as.double(weights), # vector now
            changes = as.double(changes), # matrix with n columns
            unitdistances = as.double(unitdistances), # matrix with n columns
            maxdists = as.double(maxdists),
            nobjects = as.integer(nobjects),
            nmat = as.integer(nmat), # number of data matrices
            nvar = as.integer(nvar), # vector!
            nNA = as.integer(nNA), # matrix containing the number of NAs
            ncodes = as.integer(ng),
            rlen = as.integer(rlen),
            NAOK = TRUE,
            PACKAGE = "kohonen")

  changes <- matrix(res$changes, ncol=nmat)
  colnames(changes) <- names(data)
  
  codes2 <- matrix(res$codes, nrow=ng)
  
  codes <- vector("list", length(orig.data))
  endings <- cumsum(nvar)
  beginnings <- c(1, cumsum(nvar)[-nmat] + 1)

  for (i in 1:nmat) {
    codes[[ whatmap[i] ]] <- codes2[, beginnings[i]:endings[i], drop=FALSE]
    
    if (is.factor(orig.data[[ whatmap[i] ]])) {
      colnames(codes[[ whatmap[i] ]]) <- levels(orig.data[[ whatmap[i] ]])
    } else {
      colnames(codes[[ whatmap[i] ]]) <- colnames(data[[i]])
    } 
  }
  
  names(codes) <- names(orig.data)

  if (keep.data) {
    mapping <- map.kohonen(list(codes = codes,
                                weights = orig.weights,
                                whatmap = whatmap), 
                           newdata = orig.data)
        
    structure(list(data = orig.data,
                   contin = contin,
                   na.rows = narows,
                   unit.classif = mapping$unit.classif,
                   distances = mapping$distances,
                   grid = grid,
                   codes = codes,
                   changes = changes,
                   alpha = alpha,
                   radius = radius,
                   toroidal = toroidal,
                   weights = orig.weights,
                   whatmap = whatmap,
                   method = "supersom"),
              class = "kohonen")
  } else {
    structure(list(grid = grid,
                   contin = contin,
                   codes = codes,
                   changes = changes,
                   alpha = alpha,
                   radius = radius,
                   toroidal = toroidal,
                   weights = orig.weights, 
                   whatmap = whatmap,
                   method = "supersom"),
              class = "kohonen")
  }
}

