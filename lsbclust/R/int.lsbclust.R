#' Interaction Clustering in Least Squares Bilinear Clustering
#' 
#' This function implements the interaction clustering part of the Least Squares Bilinear Clustering
#' method of Schoonees, Groenen and Van de Velden (2014). 
#' 
#' @param data A three-way array representing the data.
#' @param margin An integer giving the single subscript of \code{data} over which the clustering 
#' will be applied. 
#' @param delta A four-element binary vector (logical or numeric) indicating which sum-to-zero 
#' constraints must be enforced.
#' @param nclust An integer giving the desired number of clusters. If it is a vector, the algorithm
#' will be run for each element.
#' @param ndim The required rank for the approximation of the interactions (a scalar).
#' @param fixed One of \code{"none"}, \code{"rows"} or \code{"columns"} indicating whether to fix neither
#' sets of coordinates, or whether to fix the row or column coordinates across clusters respectively.
#' If a vector is supplied, only the first element will be used.
#' @param nstart The number of random starts to use.
#' @param starts A list containing starting configurations for the cluster membership vector. If not
#' supplied, random initializations will be generated.
#' @param alpha Numeric value in [0, 1] which determines how the singular values are distributed
#' between rows and columns.
#' @param parallelize Logical indicating whether to parallelize over different starts or not.
#' @param maxit The maximum number of iterations allowed.
#' @param verbose Integer controlling the amount of information printed: 0 = no information, 
#' 1 = Information on random starts and progress, and 2 = information is printed after
#' each iteration for the interaction clustering.
#' @param method The method for calculating cluster agreement across random starts, passed on
#' to \code{\link{cl_agreement}}. None is calculated when set to \code{NULL}.
#' @return An object of class \code{int.lsb}
#' @import parallel
#' @export
#' @examples
#' data("supermarkets")
#' out <- int.lsbclust(data = supermarkets, margin = 3, delta = c(1,1,0,0), nclust = 4, ndim = 2, 
#'            fixed = "rows", nstart = 1, alpha = 0)
#' @export
int.lsbclust <- function(data, margin = 3L, delta, nclust, ndim = 2, fixed = c("none", "rows", "columns"), 
                         nstart = 50, starts = NULL, alpha = 0.5, parallelize = FALSE, maxit = 100, 
                         verbose = 1, method = "diag"){
  
  ## Capture call, start time, original data
  time0 <- proc.time()[3]
  cll <- match.call()
  data.org <- data
  
  ## Recurse if nclust is a vector
  if (length(nclust) > 1){
    return(lapply(nclust, int.lsbclust, data = data, margin = margin, delta = delta, 
                  ndim = ndim, fixed = fixed, nstart = nstart, starts = starts, alpha = alpha,
                  parallelize = parallelize, maxit = maxit, verbose = verbose, method = method))
  }
  
  ## Sanity checks and coercion
  if (!is(data, "array") || length(dim(data)) != 3) stop("Data must be a three-way array.") 
  if (!all(margin %in% 1:3) || length(margin) != 1) stop("Argument 'margin' must be 1, 2 or 3.")
  delta <- as.numeric(delta)
  if (alpha < 0 || alpha > 1) stop("Alpha must be between 0 and 1.")
  if (length(delta) != 4 || !all(delta %in% 0:1))  stop("Argument 'delta' supplied in an unrecognized format.")
  fixed <- match.arg(tolower(fixed[1L]), choices = c("none", "rows", "columns"))
  
  ## Permute the dimensions of the data so that clustering is over the last dimension
  data <- aperm.default(data, perm = c((1:3)[-margin], margin))
  
  ## Get dimensions
  dims <- dim(data)
  J <- dims[1]
  K <- dims[2]
  N <- dims[3]
  
  ## Handle nclust == 1
  if (nclust == 1) nstart <- 1
  
  ## Check for reasonable rank
  if (ndim > min(J, K)) stop("Number of dimensions required exceeds min(J, K).")

  ## Generate random starts for cluster membership if not supplied
  if (is.null(starts)){
    starts <- lapply(X = replicate(nstart, runif(nclust) + 5/N, simplify = FALSE), 
                     FUN = sample.int, n = nclust, size = N, replace = TRUE)
  } else {
    starts <- as.list(starts)
    nstart <- length(starts)
    if(any(sapply(starts, length) != N)) 
          stop("At least one of the supplied starting configurations are of an incorrect length.")
  }
  
  ## Centre the data if necessary
  if (delta[1] == 1) {
    data <- sweep(data, MARGIN = c(2, 3), STATS = colMeans(data), FUN = "-")
  } 
  if (delta[2] == 1) {
    data <- sweep(data, MARGIN = c(1, 3), STATS = apply(data, c(1, 3), mean), FUN = "-")
  }

  ## Concatenated data
  datmat <- apply(data, 3, c)
    
  ## Function to compute loss
#   lossfun <- function(y, z) sum((y - z)^2)
  
  ## Inline function to do update of row and column coordinates given cluster membership
  updCD <- function(start, fixed){
    
    ## Cluster sizes
    nvec <- tabulate(start)
    
#     ## Calculate cluster means and cluster means weighted by the square root of cluster sizes
#     indlst <- split(seq_len(N), f = start)
#     mns <- lapply(indlst, function(x) apply(data[, , x], MARGIN = 1:2, FUN = mean))
#     ## TODO: This can be optimized
#     
#     ## Add weights if fixed != "none"
#     if(fixed != "none") mns <- mapply("*", mns, sqrt(nvec), SIMPLIFY = FALSE)
    
    ## Calculate cluster means and cluster means weighted by the square root of cluster sizes
    mns <- ClustMeans(nclust, start, datmat)
    if(fixed != "none") mns <- mns * sqrt(nvec)
    mns <- lapply(seq_len(nclust), function(x) matrix(mns[x, ], nrow = J, ncol = K))
    
    ## Function for fixed = "rows"
    CDu <- function(){
      
      ## Setup matrix concatenating along the columns
      Xc <- do.call(cbind, mns)
      
      ## Do SVD 
      svdX <- svd(Xc, nu = ndim, nv = ndim)
      
      ## Determine updates from SVD
      C <- svdX$u %*% diag(svdX$d[1:ndim]^alpha)
      Dstar <- svdX$v %*% diag(svdX$d[1:ndim]^(1 - alpha))
      
      ## Rescale Dstar and divide up in Du's
      Dstar.rs <- Dstar * 1/sqrt(nvec) %x% matrix(1, K, ndim)
      Ds <- split.data.frame(Dstar.rs, rep(1:nclust, each = K))
      
      ## Return C as Ds
      out <- list(C = list(C), D = Ds, svd = list(svdX), Xc = Xc, Xr = NULL,
                  Cstar = NULL, Dstar = Dstar, means = mns)
      return(out)
    }
      
    ## Function for fixed = "columns"
    CuD <- function(){
      
      ## Setup matrix concatenating along the rows
      Xr <- do.call(rbind, mns)
      
      ## Do SVD
      svdX <- svd(Xr, nu = ndim, nv = ndim)
      
      ## Determine updates from SVD
      Cstar <- svdX$u %*% diag(svdX$d[1:ndim]^alpha)
      D <- svdX$v %*% diag(svdX$d[1:ndim]^(1 - alpha))
      
      ## Rescale Cstar and divide up in Cu's
      Cstar.rs <- Cstar * 1/sqrt(nvec) %x% matrix(1, J, ndim)
      Cs <- split.data.frame(Cstar.rs, rep(1:nclust, each = J))
      
      ## Return C as Ds
      out <- list(C = Cs, D = list(D), svd = list(svdX), Xc =  NULL, Xr = Xr,
                  Cstar = Cstar, Dstar = NULL, means = mns)
      return(out)
    }
      
    ## Function for fixed = "none"
    CuDu <- function(){
      
      ## SVD's for all u
      svdX <- lapply(mns, svd, nu = ndim, nv = ndim)
      
      ## Calculate Cu's and Du's
      Cs <- lapply(svdX, function(x) x$u %*% diag(x$d[1:ndim]^alpha))
      Ds <- lapply(svdX, function(x) x$v %*% diag(x$d[1:ndim]^(1 - alpha)))
      
      ## Return updates of Cs and Ds
      out <- list(C = Cs, D = Ds, svd = svdX, Xc = NULL, Xr = NULL, 
                  Cstar = NULL, Dstar = NULL, means = mns)
      return(out)      
    }
    
    updfun <- switch(fixed, none = CuDu, rows = CDu, columns = CuD)
    updfun()
  }
  
  ## Function to update G, using list of model mean estimates
  updG <- function(mns){
    
    ## Matrix of loss function values
    mnsmat <- vapply(mns, c, numeric(J * K))
    lossmat <- LossMat(datmat, mnsmat)
#     lossmat <- apply(data, 3, function(x) vapply(mns, lossfun, numeric(1), z = x))
#     if(nclust == 1) lossmat <- matrix(lossmat, nrow = 1, ncol = N)
#     ## Construct arrays of the means
#     mns.arr <- lapply(mns, replicate, n = N)
#     lossmat <- mapply(FUN = lossfun, list(data), mns.arr)
    
    ## Determine class of minimum loss and the value of the loss currently
    newclass <- apply(lossmat, 2, which.min)
    losscomps <- apply(lossmat, 2, min)
    
    ## Function to check for empty classes, and if any to ensure at least on row per class
    checkempty <- function(class) {
      
        ## Class counts and determine which are empty
        classcts <- table(factor(class, levels = 1:nclust))
        if(any(classcts == 0)) {
          zeroclass <- (1:nclust)[classcts == 0]
          
          ## Move one row per empty class
          mvs <- sapply(zeroclass, function(x) which.min(lossmat[x, ]))
          class[mvs] <- zeroclass
          
          ## Update loss components and recalculate counts
          losscomps[mvs] <- lossmat[cbind(zeroclass, mvs)]
          classcts <- table(factor(class, levels = 1:nclust))
          
          ## Give message about cluster
          message("Note: an empty cluster was re-initialized.\n")
            
          ## If an empty class still exists, fill in with random choices
          if(any(classcts == 0)) {
            zeroclass <- (1:nclust)[classcts == 0]
            mvs <- sample.int(n = N, size = length(zeroclass), replace = FALSE)
            class[mvs] <- zeroclass
            losscomps[mvs] <- lossmat[cbind(zeroclass, mvs)]
            message("Note: an empty cluster was randomly re-initialized.\n")
          } 
        }
        return(class)
      }
    
    ## Update classification if empty classes and calculate loss
    newclass <- checkempty(class = newclass)
    loss <- sum(losscomps)

    return(list(newclass = newclass, loss = loss, losscomps = losscomps))
  }
  
  ## Maximum loss for standardization
  maxloss <- sum(data^2)
  
  ## Function for algorithm: loop over updCD and updG until convergence for a single start
  onestart <- function(start){
    
    ## Terminate start if one or more classes are empty    
    if(any(table(factor(start, levels = 1:nclust)) == 0)) {
      message("A random start was discarded due to empty clusters.\n")
      return(list(minloss = NA))
    }
    
    ## Initialization
    iter <- 0 
    loss <- rep(NA, maxit)
    nchange <- rep(NA, maxit)
    curclass <- start
    
    ## Loop over updCD and updG until maxit is reached or no more moves are made
    repeat {
      
      ## Increase iterations
      iter <- iter + 1
      
      ## Update Cs and D
      CDs <- updCD(start = curclass, fixed = fixed)
      
      ## Model estimates of means
      rmns <- mapply(FUN = tcrossprod, CDs$C, CDs$D, SIMPLIFY = FALSE)
      
      ## Update G
      Gupd <- updG(rmns)
      
      ## Calculate number of changes and update curclass
      ord <- clue::solve_LSAP(table(Gupd$newclass, curclass), maximum = TRUE)
      Gupd$newclass <- plyr::mapvalues(Gupd$newclass, from = 1:nclust, to = ord)
      nchange[iter] <- sum(Gupd$newclass != curclass)
      curclass <- Gupd$newclass
      
      ## Calculate standardized loss and print information if relevant
      loss[iter] <- Gupd$loss/maxloss
      if(verbose > 1) cat(sprintf(paste0("%", nchar(as.character(maxit)), "d"), iter), 
                      "| Loss =", sprintf("%5f", loss[iter]), 
                      sprintf(paste0("%", nchar(as.character(N)) + 3, "s"), 
                              paste0("(", nchange[iter], ")")),
                      "\n")
      
      ## Break in case no further changes occur or if maxit is reached
      if(iter == maxit || nchange[iter] == 0) {
        if(verbose == 1) cat("\nLoss =", sprintf("%6f", loss[iter]))
        break
      }
    }
    res <- list(C = CDs$C, D = CDs$D, svd = CDs$svd, cluster = curclass, loss = loss[1:iter], 
                nchange = nchange[1:iter], iter = iter, minloss = loss[iter], means = CDs$means, 
                Xr = CDs$Xr, Xc = CDs$Xc, Cstar = CDs$Cstar, Dstar = CDs$Dstar, 
                losscomps = Gupd$losscomps)
    class(res) <- c("int.lsbclust", "list")
#     gc(verbose = FALSE)
    return(res)
  }

  ## Apply algorithm over all starts, possibly in parallel
  if (parallelize) {
    if (.Platform$OS.type == "windows") {
      cl <- makeCluster(detectCores())
      res <- parLapply(cl = cl, X = starts, fun = onestart)
      stopCluster(cl)
    } else {
      res <- mclapply(X = starts, FUN = onestart, mc.cores = detectCores(), mc.allow.recursive = FALSE)
    }
  } else {
    res <- lapply(starts, onestart)
  }
  
  ## Determine loss for all starts
  allloss <- sapply(res, "[[", "minloss")
  
  ## Drop discarded starts
  res <- res[!is.na(allloss)]
  
  ## Determine best start
  allloss <- allloss[!is.na(allloss)]
  bestres <- res[[which.min(allloss)]]
  
  ## Calculate the number of estimated parameters 
  df <- switch(fixed, none = nclust * ndim * (J + K - ndim - sum(delta[1:2])) + N * (nclust - 1), 
               rows = ndim * (J + nclust * K - ndim - delta[1] - nclust * delta[2]) + N * (nclust - 1), 
               columns = ndim * (nclust * J + K - ndim - nclust * delta[1] - delta[2]) + N * (nclust - 1))
  
  ## Reorder class according to size
  ord <- order(table(bestres$cluster), decreasing = TRUE)
  bestres$cluster <- plyr::mapvalues(x = bestres$cluster, from = ord, to = seq_along(ord))
  if (fixed == "rows" || fixed == "none") {
    bestres$D <- bestres$D[ord]
    names(bestres$D) <- seq_len(nclust)
  }
  if (fixed == "columns" || fixed == "none") {
    bestres$C <- bestres$C[ord]
    names(bestres$C) <- seq_len(nclust)
  }
  
  ## Reorder svds
  if (fixed == "none") {
    bestres$svd <- bestres$svd[ord]
    names(bestres$svd) <- seq_len(nclust)
  }
  if (fixed == "rows") {
    inds <- unlist(split(seq_len(nclust * K), rep(seq_len(nclust), each = K))[ord], 
                   use.names = FALSE, recursive = FALSE)
    bestres$svd[[1]]$v <- bestres$svd[[1]]$v[inds, , drop = FALSE]
  }
  if (fixed == "columns") {
    inds <- unlist(split(seq_len(nclust * J), rep(seq_len(nclust), each = J))[ord], 
                   use.names = FALSE, recursive = FALSE)
    bestres$svd[[1]]$u <- bestres$svd[[1]]$u[inds, , drop = FALSE]
  }
  
  ## Reorder means, Xr, Xc, Cstar, Dstar (or a subset of these)
  bestres$means <- bestres$means[ord]
  if (fixed == "columns") {
    bestres$Xr <- do.call(rbind, split.data.frame(bestres$Xr, f = rep(seq_len(nclust), each = J))[ord])
    bestres$Cstar <- do.call(rbind, split.data.frame(bestres$Cstar, f = rep(seq_len(nclust), each = J))[ord])
  }
  if (fixed == "rows") {
    bestres$Xc <- t(do.call(rbind, split.data.frame(t(bestres$Xc), f = rep(seq_len(nclust), each = K))[ord]))
    bestres$Dstar <- do.call(rbind, split.data.frame(bestres$Dstar, f = rep(seq_len(nclust), each = K))[ord])
  }
  
  ## Calculate fit measures for biplots and persons
  if (fixed == "none") {
    CxD.comp <- Map(function(x, y) lapply(seq_len(ndim), function(z) tcrossprod(x[, z], y[, z])), bestres$C, bestres$D)
    CxD <- lapply(CxD.comp, Reduce, f = "+")
    
    ## Row fits
    rnum <- vapply(CxD, function(x) diag(tcrossprod(x)), numeric(J))
    rnum.comp <- rapply(CxD.comp, function(y) diag(tcrossprod(y)), how = "replace")
    rdenom <- vapply(bestres$means, function(x) diag(tcrossprod(x)), numeric(J))
    rnum.comp <- aperm(array(unlist(rnum.comp), dim = c(J, ndim, nclust)), c(1, 3, 2))
    rfit.comp <- sweep(rnum.comp, MARGIN = 1:2, STATS = rdenom, FUN = "/")
    rfit <- rnum/rdenom
    colnames(rfit) <- seq_len(nclust)
    rownames(rfit) <- rownames(data)
    dimnames(rfit.comp) <- list(rownames(data), seq_len(nclust), paste("Dim.", seq_len(ndim)))
    
    ## Column fits
    cnum <- vapply(CxD, function(x) diag(crossprod(x)), numeric(K))
    cnum.comp <- rapply(CxD.comp, function(y) diag(crossprod(y)), how = "replace")
    cdenom <- vapply(bestres$means, function(x) diag(crossprod(x)), numeric(K))
    cnum.comp <- aperm(array(unlist(cnum.comp), dim = c(K, ndim, nclust)), c(1, 3, 2))
    cfit.comp <- sweep(cnum.comp, MARGIN = 1:2, STATS = cdenom, FUN = "/")
    cfit <- cnum/cdenom
    colnames(cfit) <- seq_len(nclust)
    rownames(cfit) <- colnames(data)
    dimnames(cfit.comp) <- list(colnames(data), seq_len(nclust), paste("Dim.", seq_len(ndim)))
    
    ## Overall fit
    ofit <- vapply(bestres$svd, function(x) x$d^2/sum(x$d^2), numeric(min(J, K)))
    rownames(ofit) <- seq_len(min(J, K))
    
  }
  if (fixed == "columns") {
    
    ## Calculate CstarD and numerator per dimension
    CstarD.comp <- lapply(seq_len(ndim), function(x) tcrossprod(bestres$Cstar[, x], bestres$D[[1]][, x]))
    CstarD <- Reduce("+", CstarD.comp)
    rnum.comp <- sapply(CstarD.comp, function(x) diag(tcrossprod(x)))
    
    ## Row fits
    rnum <- rowSums(rnum.comp)
    rdenom <- diag(tcrossprod(bestres$Xr))
    rfit <- matrix(rnum/rdenom, ncol = nclust, nrow = J)
    rfit.comp <- array(sweep(x = rnum.comp, MARGIN = 1, STATS = rdenom, FUN = "/"), dim = c(J, nclust, ndim))
    colnames(rfit) <- seq_len(nclust)
    rownames(rfit) <- rownames(data)
    dimnames(rfit.comp) <- list(rownames(data), seq_len(nclust), paste("Dim.", seq_len(ndim)))
    
    ## Column fits
    cnum.comp <- sapply(CstarD.comp, function(x) diag(crossprod(x)))
    cnum <- rowSums(cnum.comp)
    cdenom <- diag(crossprod(bestres$Xr))
    cfit <- matrix(cnum/cdenom, ncol = 1)
    cfit.comp <- sweep(x = cnum.comp, MARGIN = 1, STATS = cdenom, FUN = "/")
    colnames(cfit) <- 1
    rownames(cfit) <- rownames(cfit.comp) <- colnames(data)
    colnames(cfit.comp) <- paste("Dim.", seq_len(ndim))
    
    ## Overall fit
    ofit <- matrix(bestres$svd[[1]]$d^2 / sum(bestres$svd[[1]]$d^2), ncol = 1)
    rownames(ofit) <- seq_len(min(nclust * J, K))
    colnames(ofit) <- 1
    
  }
  if (fixed == "rows") {
    CDstar.comp <- lapply(seq_len(ndim), function(x) tcrossprod(bestres$C[[1]][, x], bestres$Dstar[, x]))
    CDstar <- Reduce("+", CDstar.comp)
    
    ## Row fits
    rnum.comp <- sapply(CDstar.comp, function(x) diag(tcrossprod(x)))
    rnum <- rowSums(rnum.comp) # diag(tcrossprod(CDstar)) #rowSums(CDstar^2)
    rdenom <- diag(tcrossprod(bestres$Xc)) #rowSums(bestres$Xc^2)
    rfit <- matrix(rnum/rdenom, ncol = 1)
    rfit.comp <- sweep(x = rnum.comp, MARGIN = 1, STATS = rdenom, FUN = "/")
    colnames(rfit) <- 1
    rownames(rfit) <- rownames(rfit.comp) <- rownames(data)
    colnames(rfit.comp) <- paste("Dim.", seq_len(ndim))
    
    ## Column fits
    cnum.comp <- sapply(CDstar.comp, function(x) diag(crossprod(x)))
    cnum <- rowSums(cnum.comp)
    cdenom <- diag(crossprod(bestres$Xc)) #colSums(bestres$Xc^2)
    cfit <- matrix(cnum/cdenom, nrow = K, ncol = nclust)
    cfit.comp <- array(sweep(cnum.comp, MARGIN = 1, STATS = cdenom, FUN = "/"), dim = c(K, nclust, ndim))
    colnames(cfit) <- seq_len(nclust)
    rownames(cfit) <- colnames(data)
    dimnames(cfit.comp) <- list(colnames(data), seq_len(nclust), paste("Dim.", seq_len(ndim)))
    
    ## Overall fit
    ofit <- matrix(bestres$svd[[1]]$d^2/sum(bestres$svd[[1]]$d^2), ncol = 1)
    rownames(ofit) <- seq_len(min(J, nclust*K))
    colnames(ofit) <- 1
    
  }

  ## Person fit (RV/congruence coefficient)
  if (fixed != "none") CxD <- Map(tcrossprod, bestres$C, bestres$D)
  ss.X <- apply(data^2, MARGIN = 3, sum)
  ss.CD <- vapply(CxD, FUN = function(x) sum(x^2), FUN.VALUE = numeric(1))[bestres$cluster]
  ss.XCD <- rep(NA, N)
  for (i in seq_len(N)) ss.XCD[i] <- sum(diag(tcrossprod(data[, , i], CxD[[bestres$cluster[i]]])))
  pfit <- ss.XCD / sqrt(ss.X * ss.CD)
  
  ## Remove sqrt(N_u) from means if added before
  if (fixed != "none") {
    nvec <- table(bestres$cluster)
    bestres$means  <- Map("/", bestres$means, sqrt(nvec))
  }
  
  ## Update result and return
  for (i in seq_len(nclust)) {
    colnames(bestres$means[[i]]) <- colnames(data)
    rownames(bestres$means[[i]]) <- rownames(data)
  }
  names(bestres$means) <- seq_len(nclust)
  bestres$C <- lapply(bestres$C, 'rownames<-', rownames(data))
  bestres$D <- lapply(bestres$D, 'rownames<-', colnames(data))
  bestres$allloss <- allloss
  bestres$alpha <- alpha
  bestres$df <- df
  bestres$data <- data.org
  bestres$delta <- delta
  bestres$fixed <- fixed
  bestres$rfit <- rfit
  bestres$rfit.comp <- rfit.comp
  bestres$cfit <- cfit
  bestres$cfit.comp <- cfit.comp
  bestres$ofit <- ofit
  bestres$pfit <- pfit
  bestres$maxloss <- maxloss
  bestres$losscomps <- bestres$losscomps/maxloss
  bestres$nclust <- nclust
  bestres$N <- N
  bestres$call <- cll
  bestres$time <- proc.time()[3] - time0
#   bestres$means <- lapply(bestres$means, 'dimnames<-', list(rownames(data), colnames(data)))
  
  ## Cluster correspondence for best results to all other starts (if more than one start)
  if(!is.null(method) && nstart > 1) bestres$cl_agreement <- cl_agreement(res[-which.min(allloss)], 
                                                                          bestres, method = method)
  
  return(bestres)
}

#' @rdname cmat
#' @title Centring Matrix
#' @description A utility function for calculating centring matrices.
#' @param k An integer determining the dimensions of the centring matrix.
cmat <- function(k) diag(k) - matrix(1, nrow = k, ncol = k)/k

# #' @rdname centre
# #' @title Centre the Rows or Columns of a Matrix
# #' @description Remove either the column (\code{pre = TRUE}) or row (\code{pre = FALSE}) means of
# #' a matrix.
# #' @param x A matrix.
# #' @param pre Logical indicating whether to remove the column means (\code{TRUE}) or the row means
# #' (\code{FALSE}).
# # centre <- function(x, pre = TRUE) {
# #   if (!pre) {
# #     return(t(scale(t(x), scale = FALSE)))
# #   } else {
# #     return(scale(x, scale = FALSE))
# #   }
# # }
# centre <- function(x, type = c("pre", "post", "both")) {
#   type <- match.arg(type, choices = c("pre", "post", "both"))
#   retval <- switch(type,
#                    pre = x - tcrossprod(rowMeans(x), rep(1, ncol(x))),
#                    post = x - tcrossprod(rep(1, nrow(x)), colMeans(x)), 
#                    both = x - tcrossprod(rowMeans(x), rep(1, ncol(x))) - 
#                      tcrossprod(rep(1, nrow(x)), colMeans(x)) + mean.default(x))
#   return(retval)
# }
# precentre <- function(x) {
#   return(x - tcrossprod(rowMeans(x), rep(1, ncol(x))))
# }
# postcentre <- function(x) {
#   return(x - tcrossprod(rep(1, nrow(x)), colMeans(x)))
# }
# bothcentre <- function(x) {
#   x - tcrossprod(rowMeans(x), rep(1, ncol(x))) - tcrossprod(rep(1, nrow(x)), colMeans(x)) + mean.default(x)
# }