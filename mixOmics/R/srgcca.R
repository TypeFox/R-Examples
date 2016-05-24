# Copyright (C) 2015
# Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, University of Queensland, Brisbane, QLD.
# Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD

# These functions were borrowed from the RGCCA package and improved for mixOmics.

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


#############################################################################################################
# Functions modified from RGCCA
#############################################################################################################


srgcca = function (blocks, indY = NULL,  design = 1 - diag(length(blocks)), tau = NULL,
                    ncomp = rep(1, length(blocks)), scheme = "centroid", scale = TRUE,  
                    bias = FALSE, init = "svd.single", tol = .Machine$double.eps, verbose = FALSE,
                    mode = "canonical", max.iter = 500, keep = NULL, near.zero.var = FALSE, 
                    penalty = NULL) { 
  
  # blocks: list of matrices
  # indY: integer, pointer to one of the matrices of A
  # design: design matrix, links between matrices. Diagonal must be 0
  # ncomp: vector of ncomp, per matrix
  # scheme: a function "g", refer to the article (thanks Benoit)
  # scale: do you want to scale ? mean is done by default and cannot be changed (so far)
  # bias: scale the data with n or n-1
  # init: one of "svd" or "random", initialisation of the algorithm
  # tol: nobody cares about this
  # verbose: show the progress of the algorithm
  # mode: canonical, classic, invariant, regression
  # max.iter: nobody cares about this
  # keep: keepX of spls for each matrix of A. must be a list. Each entry must be of the same length (max ncomp)
  # near.zero.var: do you want to remove variables with very small variance  

  #-- checking general input parameters --------------------------------------#
  #---------------------------------------------------------------------------#
  
  #-- check that the user did not enter extra arguments
  arg.call = match.call()
  user.arg = names(arg.call)[-1]
  
  err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())),  error = function(e) e)
  
  if ("simpleError" %in% class(err))
    stop(err[[1]], ".", call. = FALSE)
  
  #-- indY
  if (!is.null(indY)){
    indY = round(indY)
    if (indY > length(blocks))
      stop("Incorrect indY")
  }
  
  
  #-- blocks
  if (!is.list(blocks))
    stop("'blocks' must be a list containing the data sets.", call. = FALSE)
  
  for (i in 1 : length(blocks)) {    
    if (is.data.frame(blocks[[i]])) 
      blocks[[i]] = as.matrix(blocks[[i]])
    
    if (!is.matrix(blocks[[i]]) || any(!is.finite(blocks[[i]])))
      stop("'block ", i, "' must be a numeric matrix.", call. = FALSE)
  }
  
  #-- equal number of rows in all blocks
  num.rows = sapply(blocks, nrow)   
  
  if (any(num.rows != nrow(blocks[[1]]))) 
    stop("unequal number of rows in blocks.", call. = FALSE)

#   #-- Check row.names
#   if (sum(sapply(lapply(blocks, row.names), is.null)) == length(blocks)) {
#     for (i in 1 : length(blocks)) row.names(blocks[[i]]) = 1 : nrow(blocks[[i]])
#   } else if (any(!sapply(blocks, function(x) {sapply(blocks, function(y) {if (is.null(row.names(x)) || is.null(row.names(y))) {return(FALSE)} else {return(all(row.names(x) %in% row.names(y)))}})}))){
#        stop ("The row.names differ between blocks")
#   }
    
  #-- block names 
  if (is.null(names(blocks))){
    names(blocks) = paste("Block", 1 : length(blocks))
  } else if (any(names(blocks) == "",na.rm=TRUE) ){
    names(blocks)[which(names(blocks) == "")] = paste("Block", which(names(blocks) == ""))
  } else if (any(is.na(names(blocks)))){
  names(blocks)[which(is.na(names(blocks)))] = paste("Block", which(is.na(names(blocks))))
  }

  #-- put variable names in each block
  for (j in 1 : length(blocks)) {
    if (is.null(colnames(blocks[[j]]))) {
      colnames(blocks[[j]]) = paste(names(blocks)[j], 1:ncol(blocks[[j]]), sep = "_") 
    }
  }
  
  #-- ncomp
  if (length(ncomp) != length(blocks))
    stop("'ncomp' must be a numeric vector of length ", length(blocks), ".", 
         call. = FALSE)
  
  if (any(is.null(ncomp)) || any(ncomp < 1) || any(!is.finite(ncomp)))
    stop("invalid values for 'ncomp'.", call. = FALSE)
  
  ncomp = round(ncomp)    
  num.var = sapply(blocks, ncol)   
  
  if (any(ncomp - num.var > 0)) 
    stop("for each block, choose a number of components smaller than the number of variables.", 
         call. = FALSE)
  
  names(ncomp) = names(blocks)
  
  #-- tau
  if (!is.null(tau)){
    if (length(tau) == 1 && is.character(tau)) {
      if (tau != "optimal")
        stop("'tau' must be either 'optimal' or a numeric vector of length ", length(blocks), ".", call. = FALSE)
    } else {
      if (length(tau) != length(blocks))
        stop("'tau' must be a numeric vector of length ", length(blocks), ".", call. = FALSE)
      
      if (any(is.null(tau)) || any(!is.finite(tau)))
        stop("invalid values for 'tau'.", call. = FALSE)
      
      if (any(tau < 0) || any(tau > 1))
        stop("'tau' values must be between 0 and 1.", call. = FALSE)
      
      names(tau) = names(blocks)
    }
  }
  
  #-- keepX and penalty parameters
  if (!is.null(keep) & !is.null(penalty))
    stop("select only one feature selection.")
  
  if (is.null(keep) & is.null(penalty))
    penalty = rep(1, length(blocks))
  
  if (!is.null(penalty)){
    if (!is.vector(penalty) & length(penalty) != length(blocks))
      stop("invalid penalty")
  }
  
  keepA = NULL
  if (!is.null(keep)){
    if (!(is.vector(keep) || is.list(keep))){
      stop ("keep must be either a vector or a list")
    }
    if (!is.list(keep)) {
      keep = list(keep = keep)
    }
    if (!is.null(keep) & (length(keep) != length(blocks))) {
      stop("keep and blocks have different lengths.")
    }
    if (is.null(ncomp)){
      stop("keep is not NULL and ncomp is NULL.")
    }
    if (any(sapply(1:length(keep), function(x){length(keep[[x]]) > ncomp[x]}))) {
      stop("length of 'keep' must be lower or equal to ", paste(ncomp, collapse = ", "), ".")
    }
    if (any(sapply(1:length(keep), function(x){any(keep[[x]] > ncol(blocks[[x]]))}))) {
      stop("each component of 'keep' must be lower or equal to the number of columns for each respective 'X'.")
    }
    
    keepA = lapply(1 : length(blocks), function(x){rep(ncol(blocks[[x]]), max(sapply(keep, length), ncomp[x]))})
    keepA = lapply(1 : length(blocks), function(x){keepA[[x]][1 : length(keep[[x]])] = keep[[x]]; return(keepA[[x]])})
  }
  
  #-- scheme
  choices = c("horst", "factorial", "centroid")
  scheme = choices[pmatch(scheme, choices)]
  
  if (is.na(scheme)) 
    stop("'scheme' should be one of 'horst', 'factorial' or 'centroid'.", 
         call. = FALSE)
  
  #-- scale
  if (!is.logical(scale))
    stop("'scale' must be a logical constant (TRUE or FALSE).", call. = FALSE)
  
  #-- bias
  if (!is.logical(bias))
    stop("'bias' must be a logical constant (TRUE or FALSE).", call. = FALSE)
  
  #-- tol
  if (is.null(tol) || tol < 0 || !is.finite(tol))
    stop("invalid value for 'tol'.", call. = FALSE)
  
  #-- verbose
  if (!is.logical(verbose))
    stop("'verbose' must be a logical constant (TRUE or FALSE).", call. = FALSE)
  
  #-- design
if ((ncol(design) != nrow(design)) || (ncol(design) != length(blocks)))
    stop("Incorrect design")
    
  #-- end checking --#
  #------------------#

  A = blocks;
  if(near.zero.var == TRUE) {
    
    nzv.A = lapply(A,nearZeroVar)
    for (q in c(1 : length(A))[ifelse(is.null(indY), 0, -indY)]) {
      if (length(nzv.A[[q]]$Position) > 0) {
        names.remove.X=colnames(A[[q]])[nzv.A[[q]]$Position]
        A[[q]] = A[[q]][, -nzv.A[[q]]$Position,drop=FALSE]
        
        if (verbose)
          warning("Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv for problematic predictors.")
        
        if(ncol(A[[q]]) == 0) {
          stop("No more variables in X")
        }
        
        #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
        if(any(keepA[[q]]>ncol(A[[q]]))){
            ind=which(keepA[[q]]>ncol(A[[q]]))
            keepA[[q]][ind]=ncol(A[[q]])
        }
      }
    }
  }
    
  # center the data per matrix of A, scale if scale=TRUE, option bias
  A = lapply(A, function(x){scale2(x, center = TRUE, scale = scale, bias = bias)})
  ni = nrow(A[[1]]) #number of samples

  ### Start: Initialization parameters
  pjs <- sapply(A, NCOL); nb_ind <- NROW(A[[1]]);
  J <- length(A); R <- A; N <- max(ncomp) # R: residuals matrices, will be a list of length ncomp
  
  AVE_inner <- AVE_outer <- rep(NA, max(ncomp))
  defl.matrix <- AVE_X <- crit <- tau.rgcca <- list()
  P <- loadings.A <- loadings.Astar <- c <- t <- b <- variates.A <- vector("list",J)
  for (k in 1:J) t[[k]] <- variates.A[[k]] <- matrix(NA, nb_ind, N)
  for (k in 1:J) P[[k]] <- loadings.A[[k]] <- loadings.Astar[[k]]<- matrix(NA, pjs[[k]], N)
  
  defl.matrix[[1]] <- A; ndefl <- ncomp - 1;  J2 <- J-1;
  
  if (is.vector(tau)){
    tau = matrix(rep(tau, N), nrow = N, ncol = length(tau), byrow = TRUE)
  }
  ### End: Initialization parameters

  iter=NULL
  for (n in 1 : N) {
    
    if (verbose)
      cat(paste0("Computation of the SGCCA block components #", n, " is under progress... \n"))
    meta.block.result = NULL
    
    ### Start: Estimation of loadings vectors
    meta.block.result <- sparse.rgcca_iteration(R, design, 
                                                tau = if (is.matrix(tau)){tau[n, ]} else {"optimal"}, 
                                                keepA = if (!is.null(keepA)) {lapply(keepA, function(x){x[n]})} else {NULL},
                                                penalty = if (!is.null(penalty)) {penalty} else {NULL},
                                                verbose = verbose, max.iter = max.iter, scheme = scheme, init = init, tol = tol)
    ### End: Estimation of loadings vectors
    
    AVE_inner[n] <- meta.block.result$AVE_inner
    crit[[n]] <- meta.block.result$crit
    tau.rgcca[[n]] <- meta.block.result$tau
    
    for (k in 1:J) 
      variates.A[[k]][, n] <- meta.block.result$variates.A[, k]
        
    # deflation if there are more than 1 component and if we haven't reach the max number of component(N)
    if (N != 1 & n != N) {
      defla.result <- defl.select(meta.block.result$variates.A, R, ndefl, n, nbloc = J, indY = indY, mode = mode, aa = meta.block.result$loadings.A)
      R <- defla.result$resdefl
      defl.matrix[[n + 1]] <- R
    }
    
    for (k in 1 : J) {
      if (N != 1  & n != N) {
        P[[k]][, n] <- defla.result$pdefl[[k]]
      }
      loadings.A[[k]][, n] <- meta.block.result$loadings.A[[k]]
    }
    
    if (n == 1) {
      for (k in 1 : J) loadings.Astar[[k]][, n] <- meta.block.result$loadings.A[[k]]
    } else {
      for (k in 1 : J) loadings.Astar[[k]][, n] <- meta.block.result$loadings.A[[k]] - loadings.Astar[[k]][, (1 : n - 1), drop = F] %*% drop(t(loadings.A[[k]][, n]) %*% P[[k]][, 1 : (n - 1), drop = F])
    }
    iter=c(iter,meta.block.result$iter)
  }
  
  if (verbose)
    cat(paste0("Computation of the SGCCA block components #", N , " is under progress...\n"))
  
  shave.matlist <- function(mat_list, nb_cols) mapply(function(m, nbcomp) m[, 1:nbcomp, drop = FALSE], mat_list, nb_cols, SIMPLIFY = FALSE)
  shave.veclist <- function(vec_list, nb_elts) mapply(function(m, nbcomp) m[1:nbcomp], vec_list, nb_elts, SIMPLIFY = FALSE)
  
  for (k in 1:J) {
    rownames(loadings.A[[k]]) = rownames(loadings.Astar[[k]])=colnames(A[[k]])
    rownames(variates.A[[k]]) = rownames(A[[k]]) #= rownames(variates.partial.A[[k]])
    colnames(variates.A[[k]]) = colnames(loadings.A[[k]]) = paste0("comp ", 1:max(ncomp))
    AVE_X[[k]] = apply(cor(A[[k]], variates.A[[k]])^2, 2, mean)
  }
  
  outer = matrix(unlist(AVE_X), nrow = max(ncomp))
  for (j in 1 : max(ncomp)) AVE_outer[j] <- sum(pjs * outer[j, ])/sum(pjs)
  variates.A = shave.matlist(variates.A, ncomp)
  AVE_X = shave.veclist(AVE_X, ncomp)
  AVE <- list(AVE_X = AVE_X, AVE_outer = AVE_outer, AVE_inner = AVE_inner)
  
  ### Start: output
  names(loadings.A)= names(variates.A) = names(A)
  names = lapply(1:J, function(x) {colnames(A[[x]])})
  names(names) = names(A)
  names[[length(names) + 1]] = row.names(A[[1]])
  names(names)[length(names)] = "indiv"
    
  ### Start: Update names list with mixOmics package
  out <- list(blocks = A, indY = indY, ncomp = ncomp,
              keep = keep, variates = variates.A, loadings = shave.matlist(loadings.A, ncomp),
              loadings.star = shave.matlist(loadings.Astar, ncomp),
              names = list(indiv = row.names(blocks[[1]]), colnames = lapply(blocks, colnames), blocks = names(blocks)),
              tol = tol, iter=iter, nzv = if(near.zero.var) nzv.A, design = design,
              scheme = scheme,  crit = crit, AVE = AVE, defl.matrix = defl.matrix, init = init, bias = bias,
              scale = scale, tau = if(!is.null(tau)) tau.rgcca, penalty = penalty)
  ### End: Update names list with mixOmics package
  
  return(out)
}

# ----------------------------------------------------------------------------------------------------------
# sparse.rgcca_iteration: this function performs the srgcca iterations until it converges for a given component
# ----------------------------------------------------------------------------------------------------------

sparse.rgcca_iteration <- function (A, design, tau = "optimal", scheme = "centroid", scale = FALSE, max.iter = 500,
                                    verbose = FALSE, init = "svd.single", bias = FALSE, tol = .Machine$double.eps, 
                                    keepA = NULL, penalty = NULL) {
  
  ### Start: Initialisation parameters
  A <- lapply(A, as.matrix)
  J <- length(A)
  n <- NROW(A[[1]])
  pjs <- sapply(A, NCOL)
  variates.A <- matrix(0, n, J)
  penalty = penalty * sqrt(pjs)
  misdata = any(sapply(A, function(x){any(is.na(x))})) # Detection of missing data
  ### End: Initialisation parameters
  
  if (!is.numeric(tau))
    tau = sapply(A, tau.estimate)
    
  loadings.A <- alpha <- M <- Minv <- K <- list()
  which.primal <- which((n >= pjs) == 1)
  which.dual <- which((n < pjs) == 1)
  
  if (init == "svd.da") {
    ### Start: Change initialization of loadings.A
     if (misdata) {
       M = lapply(c(1:(J-1)), function(x){crossprod(replace(A[[x]], is.na(A[[x]]), 0), replace(A[[J]], is.na(A[[J]]), 0))})
     } else {
       M = lapply(c(1:(J-1)), function(x){crossprod(A[[x]], A[[J]])})
     }
     
     svd.M = lapply(M, function(x){svd(x, nu = 1, nv = 1)})
     loadings.A = lapply(c(1:(J-1)), function(x){svd.M[[x]]$u})
     loadings.A[[J]] = svd.M[[1]]$v
     
     which.primal = 1 : length(A); which.dual = NULL
   } else if (init == "svd.sgcca"){
     loadings.A <- lapply(A, function(x) return(svd(x, nu = 0, nv = 1)$v))
     which.primal = 1 : length(A); which.dual = NULL
   } else if (init == "svd.rgcca") {
    for (j in which.primal) {           
      loadings.A[[j]] <- initsvd(lapply(j, function(x) {replace(A[[x]], is.na(A[[x]]), 0)})[[1]])
    }
    for (j in which.dual) {
      alpha[[j]] <- initsvd(lapply(j, function(x) {replace(A[[x]], is.na(A[[x]]), 0)})[[1]])
      K[[j]] <- A[[j]] %*% t(A[[j]])
    }
  }
  
  N = ifelse(bias, n, n - 1)
  for (j in 1 : J){
    if (j %in% which.primal) {
      M[[j]] <- ginv(tau[j] * diag(pjs[j]) + (1 - tau[j]) * cov2(A[[j]], bias = bias))    
      loadings.A[[j]] <- drop(1/sqrt(t(loadings.A[[j]]) %*% M[[j]] %*% loadings.A[[j]])) * M[[j]] %*% loadings.A[[j]]
    }
    
    if (j %in% which.dual) {
      M[[j]] = tau[j] * diag(n) + (1 - tau[j])/(N) * K[[j]]
      Minv[[j]] = ginv(M[[j]])
      alpha[[j]] = drop(1/sqrt(t(alpha[[j]]) %*% M[[j]] %*% K[[j]] %*% alpha[[j]])) * alpha[[j]]       
      loadings.A[[j]] = t(A[[j]]) %*% alpha[[j]]
    }
    
    variates.A[, j] <- A[[j]] %*% loadings.A[[j]]
  }
  
  iter = 1
  converg = crit = numeric()
  Z = matrix(0, NROW(A[[1]]), J)
  loadings.A_old = loadings.A
  g <- function(x) switch(scheme, horst = x, factorial = x^2, centroid = abs(x))
  
  repeat {
    variates.Aold <- variates.A
    
    for (j in c(which.primal, which.dual)) {
      
      if (scheme == "horst") CbyCovq <- design[j, ]
      if (scheme == "factorial") CbyCovq <- design[j, ] * cov2(variates.A, variates.A[, j], bias = bias)
      if (scheme == "centroid") CbyCovq <- design[j, ] * sign(cov2(variates.A, variates.A[, j], bias = bias))   
      
      # Compute the inner components
      Z[, j] = rowSums(mapply("*", CbyCovq, as.data.frame(variates.A)))
      
      if (j %in% which.primal) {        
        # Computer the outer weight
        loadings.A[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% M[[j]] %*% t(A[[j]]) %*% Z[, j])) * (M[[j]] %*% t(A[[j]]) %*% Z[, j])
      }
      
      if (j %in% which.dual) {        
        # Compute the outer weight
        alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j])) * (Minv[[j]] %*% Z[, j])
        loadings.A[[j]] = t(A[[j]]) %*% alpha[[j]]      
      }
    
      # sparse using keepA / penalty
      if (!is.null(keepA) || !is.null(penalty)){
        temp.norm = norm2(loadings.A[[j]])
        if (!is.null(keepA)){
          loadings.A[[j]] = sparsity(loadings.A = loadings.A[[j]], keepA = keepA[[j]], penalty = NULL)
        } else if (!is.null(penalty)){
          loadings.A[[j]] = sparsity(loadings.A = loadings.A[[j]], keepA = NULL, penalty = penalty[j])
        }
        loadings.A[[j]] = (loadings.A[[j]]/norm2(loadings.A[[j]]))*temp.norm
      }
      
      # Update variate
      variates.A[, j] = A[[j]] %*% loadings.A[[j]]
    }
    
    crit[iter] <- sum(design * g(cov2(variates.A, bias = bias)))
    
    if (iter > max.iter)
      warning(cat("The PLS algorithm did not converge after", max.iter ,"iterations."))
    
    if (max(sapply(1:J, function(x){crossprod(loadings.A[[x]] - loadings.A_old[[x]])})) < tol | iter > max.iter)
      break
    
    loadings.A_old <- loadings.A
    iter <- iter + 1
  }
    
  if (verbose) 
    plot(crit, xlab = "iteration", ylab = "criteria")
  
  AVE_inner <- sum(design * cor(variates.A)^2/2)/(sum(design)/2)
  
  result <- list(variates.A = variates.A, loadings.A = loadings.A, crit = crit[which(crit != 0)], 
                 AVE_inner = AVE_inner, design = design, tau = tau, scheme = scheme,iter=iter, keepA=keepA)
  return(result)
}

