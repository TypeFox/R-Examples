## function: overlapping group selection based on R Package 'grpreg' 
# ------------------------------------------------------------------------------
grpregOverlap <- function(X, y, group, 
                          penalty=c("grLasso", "grMCP", "grSCAD", "gel", 
                                    "cMCP", "gLasso", "gMCP"), 
                          family=c("gaussian","binomial", "poisson"), 
                          nlambda=100, lambda, 
                          lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05},
                          alpha=1, eps=.001, max.iter=1000, dfmax=ncol(X), 
                          gmax=length(group), gamma=3, tau=1/3, 
                          group.multiplier={if (strtrim(penalty,2)=="gr") 
                             sqrt(sapply(group, length)) else rep(1, length(group))}, 
                          returnX = FALSE, returnOverlap = FALSE,
                          warn=TRUE, ...) {
  # Error checking
  if (class(X) != "matrix") {
    tmp <- try(X <- as.matrix(X), silent=TRUE)
    if (class(tmp)[1] == "try-error")  {
      stop("X must be a matrix or able to be coerced to a matrix")
    }   
  }
  if (storage.mode(X)=="integer") X <- 1.0*X
  
  incid.mat <- incidenceMatrix(X, group) # group membership incidence matrix
  over.mat <- over.temp <- Matrix(incid.mat %*% t(incid.mat)) # overlap matrix
  grp.vec <- rep(1:nrow(over.mat), times = diag(over.mat)) # group index vector
  X.latent <- expandX(X, group)
  
  diag(over.temp) <- 0
  if (all(over.temp == 0)) {
    cat("Note: There are NO overlaps between groups at all!", "\n") 
    cat("      Now conducting non-overlapping group selection ...")
  }
  
  fit <- grpreg(X = X.latent, y = y, group = grp.vec, penalty=penalty,
                family=family, nlambda=nlambda, lambda=lambda, 
                lambda.min=lambda.min, alpha=alpha, eps=eps, 
                max.iter=max.iter, dfmax=dfmax, 
                gmax=gmax, gamma=gamma, tau=tau, 
                group.multiplier=group.multiplier, warn=warn, ...)
  fit$beta.latent <- fit$beta # fit$beta from grpreg is latent beta
  fit$beta <- gamma2beta(gamma = fit$beta, incid.mat, grp.vec)
  fit$incidence.mat <- incid.mat
  fit$group <- group
  fit$grp.vec <- grp.vec # this is 'group' argument in Package 'grpreg'
  if (returnX) {
    fit$X.latent <- X.latent
  } 
  if (returnOverlap) {
    fit$overlap.mat <- over.mat
  }
  # get results, store in new class 'grpregOverlap', and inherited from 'grpreg'
  val <- structure(fit,
                   class = c('grpregOverlap', 'grpreg'))
  val
}
# -------------------------------------------------------------------------------

## function: convert latent beta coefficients (gamma's) to non-latent beta's
# -------------------------------------------------------------------------------
gamma2beta<- function(gamma, incidence.mat, grp.vec) {
  # gamma: matrix, ncol = length(lambda), nrow = # of latent vars.
  p <- ncol(incidence.mat)
  J <- nrow(incidence.mat)
  beta <- matrix(0, ncol = ncol(gamma), nrow = p)
  
  intercept <- gamma[1, , drop = FALSE]
  gamma <- gamma[-1, , drop = FALSE]
  
  for (i in 1:J) {
    ind <- which(incidence.mat[i, ] == 1)
    beta[ind, ] <- beta[ind, ] + gamma[which(grp.vec == i), , drop = FALSE]
  }
  beta <- rbind(intercept, beta)         
  rownames(beta) <- c("(Intercept)", colnames(incidence.mat))
  beta
}
# -------------------------------------------------------------------------------


## function: expand design matrix X to overlapping design matrix (X.latent)
# -------------------------------------------------------------------------------
expandX <- function(X, group) {
  incidence.mat <- incidenceMatrix(X, group) # group membership incidence matrix
  over.mat <- Matrix(incidence.mat %*% t(incidence.mat), sparse = TRUE, 
                     dimnames = dimnames(incidence.mat)) # overlap matrix
  grp.vec <- rep(1:nrow(over.mat), times = diag(over.mat)) # group index vector
  
  # expand X to X.latent
  X.latent <- NULL
  names <- NULL
  
  ## TODO:
  ## (1) handle cases where variables not belong to any of groups in 'group'
  ## put each variable into a separate group, stack those group at right
  ## (2) provide option of removing groups including only one variable.
  ## Will add this later...
  
  ## the following code will automatically remove variables not included in 'group'
  for(i in 1:nrow(incidence.mat)) {
    idx <- incidence.mat[i,]==1
    X.latent <- cbind(X.latent, X[, idx, drop=FALSE])
    names <- c(names, colnames(incidence.mat)[idx])
#     colnames(X.latent) <- c(colnames(X.latent), colnames(X)[incidence.mat[i,]==1])
  }
  colnames(X.latent) <- paste('grp', grp.vec, '_', names, sep = "")
  X.latent
}
# -------------------------------------------------------------------------------
