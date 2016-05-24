# This file contains all TNAM model terms (in alphabetical order).


# outcome variable of other actors * their similarity on another attribute
attribsim <- function(y, attribute, match = FALSE, lag = 0, 
    normalization = c("no", "row", "column"), center = FALSE, coefname = NULL) {
  
  # get data in the right shape
  objects <- checkDataTypes(y = y, networks = NULL, lag = lag)
  
  # check validity of 'attribute' argument and convert
  attrib <- checkDataTypes(y = attribute, networks = NULL, lag = lag)
  if (objects$time.steps != attrib$time.steps) {
    stop("'y' and 'attribute' must have the same number of time steps.")
  }
  if (!identical(attrib$n, objects$n)) {
    stop("'attribute' must have the same number of observations as 'y'.")
  }
  attrib <- attrib$y
  
  # do the computations
  results <- list()
  for (i in 1:objects$time.steps) {
    mat <- matrix(NA, nrow = length(attrib[[i]]), ncol = length(attrib[[i]]))
    for (j in 1:length(attrib[[i]])) {
      for (k in 1:length(attrib[[i]])) {
        if (match == TRUE) {  # node match on the attribute: 1 if both the same
          if (!is.na(attrib[[i]][j]) && !is.na(attrib[[i]][k]) && 
              attrib[[i]][j] == attrib[[i]][k]) {
            mat[j, k] <- 1
          } else {
            mat[j, k] <- 0
          }
        } else if (is.na(attrib[[i]][j]) || is.na(attrib[[i]][k])) {
          mat[j, k] <- NA
        } else {  # absolute dissimilarity
          mat[j, k] <- abs(attrib[[i]][j] - attrib[[i]][k])
        }
      }  # create matrix with absolute differences on the 'attribute'
    }
    mat <- mat / max(mat, na.rm = TRUE)  # standardize
    mat <- 1 - mat  # convert to similarities
    
    # normalization of the similarity matrix
    if (normalization[1] == "row") {
      for (j in 1:nrow(mat)) {
        rs <- rowSums(mat)[j]
        for (k in 1:ncol(mat)) {
          normalized <- mat[j, k] / rs
          if (is.nan(normalized)) {
            normalized <- 0
          }
          mat[j, k] <- normalized
        }
      }
    } else if (normalization[1] == "column") {
      for (j in 1:nrow(mat)) {
        for (k in 1:ncol(mat)) {
          cs <- colSums(mat)[k]
          normalized <- mat[j, k] / cs
          if (is.nan(normalized)) {
            normalized <- 0
          }
          mat[j, k] <- normalized
        }
      }
    }
    
    results[[i]] <- mat %*% objects$y[[i]]  # apply the weights to 'y' vector
  }
  
  # convert list of results into data frame and add time and response columns
  response <- numeric()
  time <- numeric()
  attribsim <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    response <- c(response, objects$y[[i]])
    if (center == TRUE) {
      results[[i]] <- results[[i]] - mean(results[[i]], na.rm = TRUE)
    }
    attribsim <- c(attribsim, results[[i]])
  }
  
  # aggregate label
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  if (lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  label <- paste0("attribsim", coeflabel, laglabel)
  
  # aggregate and return data
  dat <- data.frame(attribsim, time = time, node = objects$nodelabels, 
      response = response)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# model term which indicates whether an actor has a certain degree centrality
centrality <- function(networks, type = c("indegree", "outdegree", "freeman", 
    "betweenness", "flow", "closeness", "eigenvector", "information", "load", 
    "bonpow"), directed = TRUE, lag = 0, rescale = FALSE, center = FALSE, 
    coefname = NULL, ...) {
  
  # check validity of arguments and prepare data
  if (is.null(directed) || !is.logical(directed)) {
    stop("'directed' must be TRUE or FALSE.")
  } else if (length(directed) != 1) {
    stop("The 'directed' argument must contain a single logical value only.")
  } else if (directed == FALSE) {
    gmode <- "graph"
  } else {
    gmode <- "digraph"
  }
  objects <- checkDataTypes(y = NULL, networks = networks, lag = lag)
  
  # do the computations
  centlist <- list()
  for (i in 1:objects$time.steps) {
    if (type[1] == "indegree") {
      cent <- degree(objects$networks[[i]], gmode = gmode, cmode = "indegree", 
          rescale = rescale, ...)
    } else if (type[1] == "outdegree") {
      cent <- degree(objects$networks[[i]], gmode = gmode, cmode = "outdegree", 
          rescale = rescale, ...)
    } else if (type[1] == "freeman") {
      cent <- degree(objects$networks[[i]], gmode = gmode, cmode = "freeman", 
          rescale = rescale, ...)
    } else if (type[1] == "betweenness") {
      cent <- betweenness(objects$networks[[i]], gmode = gmode, 
          rescale = rescale, ...)
    } else if (type[1] == "flow") {
      cent <- flowbet(objects$networks[[i]], gmode = gmode, rescale = rescale, 
          ...)
    } else if (type[1] == "closeness") {
      cent <- closeness(objects$networks[[i]], gmode = gmode, 
          rescale = rescale, ...)
    } else if (type[1] == "eigenvector") {
      cent <- evcent(objects$networks[[i]], gmode = gmode, rescale = rescale, 
          ...)
    } else if (type[1] == "information") {
      cent <- infocent(objects$networks[[i]], gmode = gmode, 
          rescale = rescale, ...)
    } else if (type[1] == "load") {
      cent <- loadcent(objects$networks[[i]], gmode = gmode, 
          rescale = rescale, ...)
    } else if (type[1] == "bonpow") {
      cent <- bonpow(objects$networks[[i]], gmode = gmode, 
          rescale = rescale, tol = 1e-20, ...)
    } else {
      stop("'type' argument was not recognized.")
    }
    centlist[[i]] <- cent
  }
  
  # aggregate data frame
  time <- numeric()
  y <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    if (is.null(centlist[[i]])) {
      y <- c(y, rep(NA, objects$n[[i]]))
    } else {
      if (center == TRUE) {
        centlist[[i]] <- centlist[[i]] - mean(centlist[[i]], na.rm = TRUE)
      }
      y <- c(y, centlist[[i]])
    }
  }
  
  # aggregate label and results and return them
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  if (lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  label <- paste0(type[1], coeflabel, laglabel)
  dat <- data.frame(y, time = time, node = objects$nodelabels)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# spatial lag for k-clique co-members
cliquelag <- function(y, networks, k.min = 2, k.max = Inf, directed = TRUE, 
    lag = 0, normalization = c("no", "row", "column"), center = FALSE, 
    coefname = NULL) {
  
  # check arguments
  if (!is.numeric(lag)) {
    stop("The 'lag' argument must be an integer value or vector of integers.")
  }
  if (is.null(k.min) || !is.numeric(k.min) || length(k.min) > 1) {
    stop("The 'k.min' argument must be a single numeric value.")
  }
  if (is.null(k.max) || !is.numeric(k.max) || length(k.max) > 1) {
    stop("The 'k.max' argument must be a single numeric value.")
  }
  if (k.max < k.min) {
    stop("'k.min' must be smaller than 'k.max'.")
  }
  if (floor(k.min) != k.min || floor(k.max) != k.max) {
    stop("The 'k.min' and 'k.max' arguments should be integer values.")
  }
  if (directed == TRUE) {
    mode <- "digraph"
  } else {
    mode <- "graph"
  }
  
  # get data in the right shape
  objects <- checkDataTypes(y = y, networks = networks, lag = lag)
  
  # do the computations
  results <- list()  # will contain vectors of results for each time step
  for (k in 1:objects$time.steps) {
    objects$network[[k]][is.na(objects$network[[k]])] <- 0  # correct NAs in nw
    objects$y[[k]][is.na(objects$y[[k]])] <- 0  # correct NAs in y
    
    # retrieve clique comembership matrices by clique size
    w3d <- clique.census(objects$networks[[k]], mode = mode, 
        clique.comembership = "bysize", tabulate.by.vertex = FALSE, 
        enumerate = FALSE)$clique.comemb
    k.max <- min(c(k.max, dim(w3d)[1]))
    w <- matrix(0, nrow = nrow(objects$networks[[k]]), 
        ncol = ncol(objects$networks[[k]]))
    for (i in k.min:k.max) {  # sum up three-cliques, four-cliques etc.
      w <- w + w3d[i, , ]
    }
    diag(w) <- 0  # diagonal is not valid (self-cliques...)
    if (normalization[1] == "row") {  # row norm. = average clique alter effect
      for (i in 1:nrow(w)) {
        rs <- rowSums(w)[i]
        for (j in 1:ncol(w)) {
          normalized <- w[i, j] / rs
          if (is.nan(normalized)) {
            normalized <- 0
          }
          w[i, j] <- normalized
        }
      }
    } else if (normalization[1] == "column") {
      for (i in 1:nrow(w)) {
        for (j in 1:ncol(w)) {
          cs <- colSums(w)[j]
          normalized <- w[i, j] / cs
          if (is.nan(normalized)) {
            normalized <- 0
          }
          w[i, j] <- normalized
        }
      }
    }
    results[[k]] <- w %*% objects$y[[k]]
  }
  
  # convert list of results into data frame and add time and response columns
  response <- numeric()
  time <- numeric()
  cliquelag <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    response <- c(response, objects$y[[i]])
    if (center == TRUE) {
      results[[i]] <- results[[i]] - mean(results[[i]], na.rm = TRUE)
    }
    cliquelag <- c(cliquelag, results[[i]])
  }
  
  # create the label
  if (length(lag) == 1 && lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  if (normalization[1] == "row") {
    normlabel <- ".rownorm"
  } else if (normalization[1] == "column" || normalization[1] == "col") {
    normlabel <- ".colnorm"
  } else {
    normlabel <- ""
  }
  klabel <- paste0(".k.", k.min, ".", k.max)
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  label <- paste0("cliquelag", coeflabel, klabel, laglabel, normlabel)
  
  dat <- data.frame(cliquelag, time = time, node = objects$nodelabels, 
      response = response)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# local clustering coefficient
clustering <- function(networks, directed = TRUE, lag = 0, center = FALSE, 
    coefname = NULL, ...) {
  
  # prepare arguments and data
  if (directed == TRUE) {
    mode <- "directed"
  } else {
    mode <- "undirected"
  }
  objects <- checkDataTypes(y = NULL, networks = networks, lag = lag)
  
  # do the computations
  results <- list()
  for (i in 1:objects$time.steps) {
    # local clustering coefficient from the igraph package
    g <- graph.adjacency(objects$networks[[i]], mode = mode)
    trans <- transitivity(g, type = "local", isolates = "zero")
    results[[i]] <- trans
  }
  
  # convert list of results into data frame and add time column
  time <- numeric()
  lcc <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    if (is.null(results[[i]])) {
      lcc <- c(lcc, rep(NA, objects$n[[i]]))
    } else {
      if (center == TRUE) {
        results[[i]] <- results[[i]] - mean(results[[i]], na.rm = TRUE)
      }
      lcc <- c(lcc, results[[i]])
    }
  }
  
  # aggregate label and return results
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  if (lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  label <- paste0("clustering", coeflabel, laglabel)
  dat <- data.frame(lcc, time = time, node = objects$nodelabels)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# simple covariate or temporally lagged covariate
covariate <- function(y, lag = 0, exponent = 1, center = FALSE, 
    coefname = NULL) {
  
  # check validity of arguments
  if (is.null(lag) || !is.numeric(lag) || floor(lag) != lag) {
    stop("'lag' should be an integer value or a vector of integers.")
  }
  if (is.null(exponent)) {
    exponent <- 1
  } else if (length(exponent) > 1 || !is.numeric(exponent)) {
    stop("'exponent' should contain a single numeric value.")
  }
  objects <- checkDataTypes(y = y, networks = NULL, lag = lag)
  
  response <- numeric()
  time <- numeric()
  y <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    response <- c(response, objects$y[[i]])
    if (center == TRUE) {
      result <- objects$y[[i]]^exponent - mean(objects$y[[i]]^exponent, 
          na.rm = TRUE)
    } else {
      result <- objects$y[[i]]^exponent
    }
    y <- c(y, result)
  }
  
  # aggregate label
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  if (length(lag) == 1 && lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  label <- paste0("covariate", coeflabel, laglabel)
  dat <- data.frame(y, time = time, node = objects$nodelabels, 
      response = response)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# model term which indicates whether an actor has a certain degree centrality
degreedummy <- function(networks, deg = 0, type = c("indegree", "outdegree", 
    "freeman"), reverse = FALSE, directed = TRUE, lag = 0, center = FALSE, 
    coefname = NULL, ...) {
  
  # check validity of arguments and prepare data
  if (is.null(directed) || !is.logical(directed)) {
    stop("'directed' must be TRUE or FALSE.")
  } else if (length(directed) != 1) {
    stop("The 'directed' argument must contain a single logical value only.")
  } else if (directed == FALSE) {
    gmode <- "graph"
  } else {
    gmode <- "digraph"
  }
  if (is.null(deg) || !is.numeric(deg) || any(deg < 0)) {
    stop("'deg' must be a numeric value or vector of positive integers.")
  }
  objects <- checkDataTypes(y = NULL, networks = networks, lag = lag)
  
  # do the computations
  dummylist <- list()
  for (i in 1:objects$time.steps) {
    d <- degree(objects$networks[[i]], gmode = gmode, cmode = type[1], ...)
    d <- 1 * (d %in% deg)
    if (reverse == TRUE) {
      d <- d * -1 + 1
    }
    dummylist[[i]] <- d
  }
  
  # aggregate data frame
  time <- numeric()
  y <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    if (is.null(dummylist[[i]])) {
      y <- c(y, rep(NA, objects$n[[i]]))
    } else {
      if (center == TRUE) {
        dummylist[[i]] <- dummylist[[i]] - mean(dummylist[[i]], na.rm = TRUE)
      }
      y <- c(y, dummylist[[i]])
    }
  }
  
  # aggregate label and return data
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  if (lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  degreelabel <- paste0(".degree", paste(deg, collapse = "."))
  label <- paste0(type[1], coeflabel, degreelabel, laglabel)
  dat <- data.frame(y, time = time, node = objects$nodelabels)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# interaction effect between multiple terms
interact <- function(x, y, lag = 0, center = FALSE, coefname = NULL) {
  
  # create interaction between x and y
  if (class(x) == "data.frame") {
    first <- x[, 1]
  } else {
    first <- x
  }
  if (class(y) == "data.frame") {
    second <- y[, 1]
  } else {
    second <- y
  }
  z <- first * second
  
  # create time and node variables
  if (class(x) == "data.frame") {
    time <- x[, 2]
    node <- as.character(x[, 3])
  } else if (class(y) == "data.frame") {
    time <- y[, 2]
    node <- as.character(y[, 3])
  } else {
    time <- rep(1, length(z))
  }
  
  # centering
  if (center == TRUE) {
    u <- unique(time)
    for (i in 1:length(u)) {
      indices <- which(time == i)
      z[indices] <- z[indices] - mean(z[indices], na.rm = TRUE)
    }
  }
  
  # create label
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  label <- paste0("interaction", coeflabel)
  
  # aggregate results, add response variable if available, and return result
  if (class(x) == "data.frame" && ncol(x) == 4) {
    dat <- data.frame(z, time, node, x[, 4])
    colnames(dat) <- c(label, "time", "node", "response")
  } else if (class(y) == "data.frame" && ncol(y) == 4) {
    dat <- data.frame(z, time, node, y[, 4])
    colnames(dat) <- c(label, "time", "node", "response")
  } else {
    dat <- data.frame(z, time, node)
    colnames(dat) <- c(label, "time", "node")
  }
  dat$node <- as.character(dat$node)
  if (is.null(lag)) {
    lag <- union(attributes(x)$lag, attributes(y)$lag)
  } else if (!is.numeric(lag[1])) {
    stop("'lag' must be a single integer value.")
  } else if (attributes(x)$lag != 0 || attributes(x)$lag != 0) {
    message(paste("Using the 'lag' argument from the interaction term, not", 
        "from the 'x' or 'y' arguments handed over to the interaction term."))
  }
  if (length(lag) > 1) {
    warning(paste("Interaction term: several 'lag' arguments are provided.", 
        "Only the first one is retained."))
  }
  if (is.null(attributes(dat)$lag)) {
    attributes(dat)$lag <- 0
  }
  return(dat)
}


# spatial network lag term with optional temporal lag and path distance decay
netlag <- function(y, networks, lag = 0, pathdist = 1, 
    decay = pathdist^-1, normalization = c("no", "row", "column", "complete"), 
    reciprocal = FALSE, center = FALSE, coefname = NULL, ...) {
  
  # check validity of arguments
  if (!is.numeric(pathdist)) {
    stop("'pathdist' must be numeric.")
  }
  if (length(pathdist) != length(decay)) {
    stop("'decay' and 'pathdist' must have the same length.")
  }
  if (any(!is.finite(decay)) || any(!is.finite(pathdist))) {
    stop("'pathdist' and/or 'decay' contain(s) infinite values.")
  }
  if (is.null(reciprocal) || length(reciprocal) > 1) {
    stop("The 'reciprocal' argument must be TRUE or FALSE.")
  }
  
  # get data in the right shape
  objects <- checkDataTypes(y = y, networks = networks, lag = lag)
  
  # do the computations
  results <- list()  # will contain vectors of results for each time step
  for (k in 1:objects$time.steps) {
    pdistmat <- geodist(objects$networks[[k]], ...)$gdist  # path dist mat.
    results[[k]] <- suppressWarnings(netLagCppLoop(objects$networks[[k]], 
        pdistmat, pathdist, decay, objects$y[[k]], normalization[1], 
        reciprocal))  # netlag loop in C++
  }
  
  # convert list of results into data frame and add time and response columns
  response <- numeric()
  time <- numeric()
  spatlag <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    response <- c(response, objects$y[[i]])
    if (center == TRUE) {
      results[[i]] <- results[[i]] - mean(results[[i]], na.rm = TRUE)
    }
    spatlag <- c(spatlag, results[[i]])
  }
  
  if (length(lag) == 1 && lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  pathdistlabel <- paste0(".pathdist", paste(pathdist, collapse = "."))
  if (any(decay != 1)) {
    decaylabel <- paste0(".decay", paste(decay, collapse = "."))
  } else {
    decaylabel <- ""
  }
  if (normalization[1] == "row") {
    normlabel <- ".rownorm"
  } else if (normalization[1] == "column" || normalization[1] == "col") {
    normlabel <- ".colnorm"
  } else {
    normlabel <- ""
  }
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  label <- paste0("netlag", coeflabel, laglabel, pathdistlabel, decaylabel, 
      normlabel)
  dat <- data.frame(spatlag, time = time, node = objects$nodelabels, 
      response = response)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# structural similarity term with optional temporal lag
structsim <- function(y, networks, lag = 0, method = c("euclidean", 
    "minkowski", "jaccard", "binary", "hamming"), center = FALSE, 
    coefname = NULL, ...) {
  
  # get data in the right shape
  objects <- checkDataTypes(y = y, networks = networks, lag = lag)
  
  # do the computations
  results <- list()  # will contain vectors of results for each time step
  for (k in 1:objects$time.steps) {
    # compute structural similarity matrix
    if (method[1] == "euclidean" || method[1] == "minkowski") {
      d <- as.matrix(dist(objects$networks[[k]], method = method[1], ...))
      mx <- max(d, na.rm = TRUE)
      s <- (mx - d) / mx  # convert dist to similarity and standardize [0; 1]
    } else if (method[1] == "binary") {
      d <- as.matrix(dist(objects$networks[[k]], method = method[1], ...))
      s <- 1 - d
    } else if (method[1] == "jaccard") {
      d <- as.matrix(vegdist(objects$networks[[k]], method = "jaccard", 
          na.rm = TRUE, ...))
      s <- 1 - d
    } else if (method[1] == "hamming") {
      d <- sedist(objects$networks[[k]], method = "hamming", ...) / 
          nrow(objects$networks[[k]])  # standardize to [0; 1]
      s <- 1 - d
    } else {
      stop("'method' argument was not recognized.")
    }
    diag(s) <- 0
    rm(d)
    
    # apply structural similarities as weight matrix
    s[is.na(s)] <- 0
    result <- s %*% objects$y[[k]]
    results[[k]] <- result
  }
  
  # convert list of results into data frame and add time and response columns
  response <- numeric()
  time <- numeric()
  structsim <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    response <- c(response, objects$y[[i]])
    if (center == TRUE) {
      results[[i]] <- results[[i]] - mean(results[[i]], na.rm = TRUE)
    }
    structsim <- c(structsim, results[[i]])
  }
  
  # aggregate label
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  if (lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  label <- paste0("structsim", coeflabel, laglabel, ".", method[1])
  
  # aggregate and return data
  dat <- data.frame(structsim, time = time, node = objects$nodelabels, 
      response = response)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# spatial weighted lag with row and column normalization; for weighted matrices
weightlag <- function(y, networks, lag = 0, normalization = c("no", "row", 
    "column"), center = FALSE, coefname = NULL) {
  
  # check arguments
  if (!is.numeric(lag)) {
    stop("The 'lag' argument must be an integer value or vector of integers.")
  }
  
  # get data in the right shape
  objects <- checkDataTypes(y = y, networks = networks, lag = lag)
  
  # do the computations
  results <- list()  # will contain vectors of results for each time step
  for (k in 1:objects$time.steps) {
    objects$networks[[k]][is.na(objects$networks[[k]])] <- 0
    if (normalization[1] == "row") {
      for (i in 1:nrow(objects$networks[[k]])) {
        rs <- rowSums(objects$networks[[k]])[i]
        for (j in 1:ncol(objects$networks[[k]])) {
          normalized <- objects$networks[[k]][i, j] / rs
          if (is.nan(normalized)) {
            normalized <- 0
          }
          objects$networks[[k]][i, j] <- normalized
        }
      }
    } else if (normalization[1] == "column") {
      for (i in 1:nrow(objects$networks[[k]])) {
        for (j in 1:ncol(objects$networks[[k]])) {
          cs <- colSums(objects$networks[[k]])[j]
          normalized <- objects$networks[[k]][i, j] / cs
          if (is.nan(normalized)) {
            normalized <- 0
          }
          objects$networks[[k]][i, j] <- normalized
        }
      }
    }
    result <- objects$networks[[k]] %*% objects$y[[k]]
    results[[k]] <- result
  }
  
  # convert list of results into data frame and add time and response columns
  response <- numeric()
  time <- numeric()
  spatlag <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    response <- c(response, objects$y[[i]])
    if (center == TRUE) {
      results[[i]] <- results[[i]] - mean(results[[i]], na.rm = TRUE)
    }
    spatlag <- c(spatlag, results[[i]])
  }
  
  # create the label
  if (length(lag) == 1 && lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  if (normalization[1] == "row") {
    normlabel <- ".rownorm"
  } else if (normalization[1] == "column" || normalization[1] == "col") {
    normlabel <- ".colnorm"
  } else {
    normlabel <- ""
  }
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  label <- paste0("weightlag", coeflabel, laglabel, normlabel)
  
  dat <- data.frame(spatlag, time = time, node = objects$nodelabels, 
      response = response)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}

