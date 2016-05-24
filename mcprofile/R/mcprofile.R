mcprofile <-
function(object, CM, control=mcprofileControl(), grid=NULL) UseMethod("mcprofile")

mcprofile.glm <-
function(object, CM, control=mcprofileControl(), grid=NULL){
  if (is.null(rownames(CM))) rownames(CM) <- paste("C",1:nrow(CM), sep="")
  if (is.null(colnames(CM))) colnames(CM) <- names(coefficients(object))
  if (ncol(CM) != length(coefficients(object))) stop("Equal number of contrast and model coefficients needed!")

  df.needed <- family(object)$family == "gaussian" | length(grep("quasi", family(object)$family)) == 1 | length(grep("Negative Binomial", family(object)$family)) == 1

  ## construct grid
  if (is.null(grid)) grid <- constructGrid(object, CM, control)
  if (is.matrix(grid)) grid <- lapply(1:ncol(grid), function(i) na.omit(grid[,i]))
  ## check grid
  if (length(grid) != nrow(CM)) stop("Number of contrasts and grid support vectors differ!")
  if (any(sapply(grid, function(g) any(rank(g) != 1:length(g))))) stop("grid support is not in increasing order!") 
  if (any(sapply(grid, length) < 2)) stop("At least 2 grid supports per contrast are needed!")  

  ## model info
  est <- coefficients(object)
  OriginalDeviance <- object$deviance
  DispersionParameter <- summary(object)$dispersion ### Dispersion estimate @ every step needed?
  mf <- model.frame(object)
  Y <- model.response(mf)
  n <- NROW(Y)
  O <- model.offset(mf)
  if (!length(O)) O <- rep(0, n)
  W <- model.weights(mf)
  if (length(W) == 0L) W <- rep(1, n)
  X <- model.matrix(object)
  fam <- family(object)
  etastart <- X %*% est
  glmcontrol <- object$control

  ## profiling
  prolist <- lapply(1:nrow(CM), function(i){
    K <- CM[i,,drop=FALSE]
    glmpro <- glm_profiling(X, Y, W, etastart, O, fam, glmcontrol, est, OriginalDeviance, DispersionParameter, K, grid[[i]])
    return(glmpro)
  })  
    
  out <- list()
  out$object <- object
  out$CM <- CM
  out$srdp <- lapply(prolist, function(x) x[[1]])
  out$optpar <- lapply(prolist, function(x) x[[2]])
  out$cobject <- lapply(prolist, function(x) x[[3]])
  if (df.needed) out$df <- df.residual(object) else df <- NULL
  class(out) <- "mcprofile"
  out
}

mcprofile.lm <-
function(object, CM, control=mcprofileControl(), grid=NULL){
  oc <- as.list(object$call)
  oc$family <- call("gaussian")
  oc[[1]] <- as.symbol("glm")
  object <- eval(as.call(oc))
  mcprofile.glm(object, CM=CM, control=control, grid=grid)
}


mcprofileControl <-
function(maxsteps=10, alpha=0.01, del=function(zmax) zmax/5){
  list(maxsteps=maxsteps, alpha=alpha, del=del)
}
