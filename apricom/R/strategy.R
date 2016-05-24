#' An Accessory Function for Strategy Selection
#'
#' This function acts as a switchboard for various regression modelling strategy,
#' and links these strategies to \code{\link{compare}}.
#'
#' For further details, see \code{\link{compare}}.
#'
#' @import logistf
#' @import penalized
#' @import rms
#' @importFrom stats glm.fit binomial update
#'
#' @param g        a data matrix to which a modelling strategy will be applied.
#' @param method    a list containing the a regression modelling strategy name
#'                and strategy-specific parameter values.
#' @param int      logical. If int is TRUE an intercept will be included in the
#'                 regression model.
#' @param int.adj  logical. If int.adj is TRUE the intercept will be re-estimated
#'                 after shrinkage is applied to the regression coefficients.
#' @param model    the type of regression model. Either "linear" or "logistic".
#'
#' @return The object returned by \code{strategy} depends on the strategy selected in the
#'         call to \code{compare}. For example, if \code{method[1] = "lsq"}, \code{strategy} will return
#'         a matrix of regression coefficients, whereas if \code{method[1] = "boot"}, strategy
#'         will return a list containing values returned by \code{bootval}.
#'
#' @note This function is not designed to be called directly, but acts as a workhorse
#'      function for \code{compare}


strategy <- function(g, method, int, int.adj, model) {

  flag <- method[1]
  if (int == TRUE) sdm0 <- matrix(c(0, rep(1, dim(g)[2] - 2)), nrow = 1)
  else sdm0 <- matrix(rep(1, dim(g)[2] - 1), nrow = 1)


### ORDINARY LEAST SQUARES ESTIMATION
  if (flag == "lsq") {
    strat.out <- ols.rgr(g)
    return(strat.out)
  }


### MAXIMUM LIKELIHOOD ESTIMATION
  if (flag == "ml") {
    strat.out <- ml.rgr(g)
    return(strat.out)
  }


### SHRINKAGE USING HEURSTIC FORMULAE
  if (flag == "heuristic") {
    if(length(method) < 2) warning("The number of degrees of freedom used by predictors was estimated automatically")
    if (length(method) < 2) DF <- dim(g)[2] - 2
     else DF <- method[[2]]
    strat.out <- shrink.heur(g, model, DF, int)
    return(strat.out)
  }


### SPLIT-SAMPLE SHRINKAGE AFTER ESTIMATION
  if (flag == "split") {
    if (length(method) < 4) {method[[4]] <- sdm0}
    if (length(method) < 3) {method[[3]] <- 0.5}
    if (length(method) < 2) {method[[2]] <- 1}
    nrounds <- method[[2]]
    fract <- method[[3]]
    sdm <- method[[4]]
    strat.out <- splitval(g, model, nrounds, fract, sdm, int, int.adj)
    return(strat.out)
  }


### BOOTSTRAP SHRINKAGE AFTER ESTIMATION
  if (flag == "boot") {
    if (length(method) < 3){method[[3]] <- sdm0}
    N <- method[[2]]
    sdm <- method[[3]]
    strat.out <- bootval(g, model, N, sdm, int, int.adj)
    return(strat.out)
  }

### CROSS-VALIDATION SHRINKAGE AFTER ESTIMATION
if (flag == "kcv") {
  if (length(method) < 4) {method[[4]] <- sdm0}
  if (length(method) < 3) {method[[3]] <- 1}
  k <- method[[2]]
  nreps <- method[[3]]
  sdm <- method[[4]]
  strat.out <- kcrossval(g, model, k, nreps, sdm, int, int.adj)
  return(strat.out)
}

### LEAVE-ONE-OUT CROSS-VALIDATION SHRINKAGE AFTER ESTIMATION
if (flag == "loocv") {
  if (length(method) < 2) {method[[2]] <- 1}
  if (length(method) < 3) {method[[3]] <- sdm0}
  nreps <- method[[2]]
  sdm <- method[[3]]

  strat.out <- loocval(g, model, nreps, sdm, int, int.adj)
  return(strat.out)
}

### PENALIZED MAXIMUM LIKELIHOOD ESTIMATION Likelihood Cross Validation
  if (flag == "pml.LCV") {

    folds = method[[2]]
    if (int == TRUE) gdf <- as.data.frame(g[, 2:dim(g)[2]])
    else gdf <- as.data.frame(g)
    gy <- gdf[, dim(gdf)[2]]
    gX <- as.matrix(gdf[, 1:dim(gdf)[2] - 1])

    if (int == TRUE) L2fit <- optL2(gy ~ gX, data = gdf, model = model, fold = folds, maxiter = 50, trace = FALSE)

    else L2fit <- optL2(gy ~ gX, unpenalized = ~ 0, data = gdf, model = model, fold = nrounds, maxiter = 50, trace = FALSE)

    strat.out <- L2fit
    return(strat.out)
  }

### PENALIZED MAXIMUM LIKELIHOOD ESTIMATION Firth's/Heinze method

  if (flag == "pml.firth") {

    if (int == TRUE) {
      g.Y <- g[, dim(g)[2]]
      g.X <- as.data.frame(g[, 3:dim(g)[2] - 1])
      y <- matrix(logistf(formula = g.Y ~ ., data = g.X)$coefficients, ncol = 1)
      return(y)
    } else {
      g.Y <- g[, dim(g)[2]]
      g.X <- as.data.frame(g[, 1:dim(g)[2]])
      strat.out <- logistf(formula = g.Y ~ . - 1, data = g.X)
      return(strat.out)
    }
  }

### PENALIZED MAXIMUM LIKELIHOOD ESTIMATION Harrell rms library, based on optimizing the corrected AIC

  if (flag == "pml.AIC") {

    if (int == FALSE) stop("For PML with a corrected AIC-derived penalty, a regression intercept must be included")


    if (length(method) < 2) penalty.range <- -100:200
    else penalty.range = method[[2]]

    g.Y <- g[, dim(g)[2]]
    g.X <- as.data.frame(g[, 3:dim(g)[2] - 1])
    Penfit <- lrm(formula = g.Y ~ ., data = g.X, x = TRUE, y = TRUE)
    PenY <- pentrace(Penfit, penalty = penalty.range)
    strat.out <- update(Penfit, penalty = PenY$penalty)
    return(strat.out)
  }

  else stop("Strategy must be of type 'lsq', 'ml', 'heuristic', 'split', 'kcv', 'boot', 'pml.LCV', 'pml.firth' or 'pml.AIC'")

}
