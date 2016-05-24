#' @include spikeSlabGAM.R
{}

#' Summary for posterior of a \code{spikeSlabGAM} fit
#'
#' Returns basic information about prior and model structure, posterior means of
#' inclusion probablities, term inportance and the most probable models found by
#' the SSVS.
#'
#' Two scalar summaries of term importance are included: \describe{
#' \item{\code{P(gamma = 1)}}{The posterior mean of \eqn{P(\gamma = 1)}, i.e. the
#' marginal posterior inclusion probability of the term.} \item{\code{pi}}{The
#' scaled dot products of the posterior mean of the linear predictor
#' \eqn{\eta_i} associated with the \eqn{i}-th term and the total linear
#' predictor \eqn{\eta}: \eqn{\pi_i = \eta_i^T \eta/ \eta^T \eta}. They sum to
#' one (but can be negative as well), and (for mutually orthogonal \eqn{\eta_i})
#' provide a percentage decomposition of the sum of squares of \eqn{\eta}. See
#' references for details.}} The summary also shows the dimensionality of the
#' basis associated with a term.
#'
#' The top row in the model table shows the relative frequency of the respective
#' model, the bottom row shows cumulative relative frequencies for the models
#' visited most often.
#' @param object an object of class \code{spikeSlabGAM}
#' @param threshold threshold for inclusion of a model term.  Defaults to 0.5.
#' @param ... 	arguments passed from or to other methods (not used)
#' @return an object of class \code{summary.spikeSlabGAM}
#' @author Fabian Scheipl
#' @references Gu, Chong (2002). Smoothing Spline ANOVA models. Springer. (see
#'   chapter 3.6)
#' @export
summary.spikeSlabGAM <- function(object, threshold =.5, ...) {
  stopifnot(class(object)=="spikeSlabGAM", threshold >= 0, threshold <= 1)

  inds <- sapply(unique(names(object$model$groupIndicatorsOrig)),
    function(x) which(names(object$model$groupIndicatorsOrig)== x),
    simplify = FALSE)

  dims <- sapply(inds, length)
  etas <- sapply(inds, function(x) object$X[, x, drop = F] %*%
      object$postMeans$beta[x, drop = F])
  q <- NCOL(etas)

  etas <- cbind(etas, eta = rowSums(etas))
  if(object$family == 0) {
    resid <- object$y - etas[,"eta"]
    etas <- cbind(etas, resid)
    colnames(etas)[NCOL(etas)] <- "resid"
  }

  trmSummary <- matrix(NA, nrow = q, ncol = 3)
  rownames(trmSummary) <- unique(names(object$model$groupIndicatorsOrig))
  colnames(trmSummary) <- c("P(gamma = 1)", "pi", "dim")
  penRows <- (rownames(trmSummary)!="u") & !grepl("u(", rownames(trmSummary),
    fixed = T)
  penCols <- if(any(!penRows)) {
    (1:q)[-which(!penRows)]
  } else 1:q

  trmSummary[penRows, "P(gamma = 1)"] <- round(object$postMeans$pV1, 3)
  etaNoU <- if(any(!penRows)) {
    etas[,"eta"] - rowSums(etas[, which(!penRows), drop = FALSE])
  } else etas[,"eta"]
  etaNoUSq <- crossprod(etaNoU)
  trmSummary[penRows, "pi"] <- apply(etas[, penCols, drop = FALSE], 2,
    function(x) t(etaNoU) %*% x/etaNoUSq)
  trmSummary[1:q, "dim"] <- dims

  response <- if(object$family == 1 & any(object$model$scale!= 1)) {
    cbind(object$y * object$model$scale, object$model$scale -
        object$y * object$model$scale)
  } else object$y
  nullDev <- -2 * logLik(glm(response ~ 1, offset = object$model$offset,
    family = switch(as.character(object$family),
        "0"= "gaussian",
        "1"= "binomial",
        "2"= "poisson")))


  modelTable <- getModels(object, threshold)


  return(structure(c(list(modelTable = modelTable, trmSummary = trmSummary,
    nullDev = nullDev, thresh = threshold), object),
    class ="summary.spikeSlabGAM"))
}

#' Print summary for posterior of a \code{spikeSlabGAM} fit
#'
#' The model table shows at least the 10 and at most the 20 most probable models
#' as found by SSVS, or enough models to account for 90\% of the probability
#' mass in the model space as found by SSVS by default.
#'
#' @param x an object of class \code{summary.spikeSlabGAM}
#' @param printPGamma  print marginal inclusion probabilities and norm of the
#'   linear predictor for each term? Defaults to TRUE
#' @param printModels  print most probable models visited by the SSVS? Defaults
#'   to TRUE
#' @param showModels how many of th most probable models to display. See details
#'   for default.
#' @param digits see \code{\link[base]{options}}, defaults to 3
#' @param ... 	arguments based from or to other methods (not used)
#' @seealso \code{\link{summary.spikeSlabGAM}}
#' @return invisibly returns \code{x}
#' @author Fabian Scheipl
#' @method print summary.spikeSlabGAM
#' @export
print.summary.spikeSlabGAM <- function(x, digits = 3, printPGamma = TRUE,
  printModels = TRUE, showModels = NULL, ...) {
  cat(paste("Spike-and-Slab STAR for", switch(as.character(x$family),
    "0"= "Gaussian",
    "1"= "Binomial",
    "2"= "Poisson"), "data \n\nModel:\n"))

  if(version$minor >= 12) {
    print(x$formula, showEnv = FALSE)
  }else{
    print(x$formula)
  }

  wrapCat(paste(x$model$n, "observations;", x$model$q, "coefficients in",
    length(unique(names(x$model$groupIndicators))), "model terms."))
  cat("\n\nPrior:\n")
  priorVec <- with(x$hyperparameter,
    c("a[tau]"= unname(tau2[1]), "b[tau]"= unname(tau2[2]),
      "v[0]"= unname(gamma["v0"]), "a[w]"= unname(w[1]), "b[w]"= unname(w[2])))
  if(x$family == 0) priorVec <- c(priorVec, with(x$hyperparameter,
    c("a[sigma^2]"= unname(sigma2[1]), "b[sigma^2]"= unname(sigma2[2]))))
  print(priorVec, digits = digits)
  cat("\nMCMC:\n")
  mcmcVec <- c(x$mcmc[1:4])
  with(x$mcmc, wrapCat(paste("Saved",  nChains * chainLength, "samples from",
    nChains, "chain(s), each ran", chainLength * thin,
    "iterations after a burn-in of", burnin, "; Thinning:", thin,"\n")))
  cat("\n")
  mcmcOpts <- with(x$mcmc, paste(ifelse(!useRandomStart,
    "Did not use randomized starting values.", ""),
    ifelse(!allKsi,
      "Didn't use parameter expansion for terms with 1 coefficient.",  ""),
    ifelse(x$family!= 0,
      paste("\nP-IWLS acceptance rates: ", round(accept[1], 2)," for alpha; ",
        round(accept[2], 2)," for xi.", sep =""),  ""),
    ifelse(x$family!= 0 && modeSwitching != 0.05,
      paste("Set 'Mode Switching'-rate to ", modeSwitching,".", sep =""),  ""),
    "\n"))

  wrapCat(mcmcOpts)
  cat("\n")

  devmat <-as.table(matrix(c("Null deviance:"= x$nullDev,
    "Mean posterior deviance:"= x$DIC["Dbar"]), nrow = 2, ncol = 1))
  rownames(devmat) <- c("Null deviance:", "Mean posterior deviance:")
  colnames(devmat) <- c("")
  print(as.table(devmat), digits = digits)
  cat("\n")

  if(printPGamma) {
    wrapCat("\nMarginal posterior inclusion probabilities and term importance:")
    cat("\n")
    stars <- paste(ifelse(x$trmSummary[,"P(gamma = 1)"]>.25,"*",""),
      ifelse(x$trmSummary[,"P(gamma = 1)"]>.5,"*",""),
      ifelse(x$trmSummary[,"P(gamma = 1)"]>.9,"*",""), sep ="")
    stars[grepl("NA", stars)] <- ""
    printTrms <- data.frame(round(x$trmSummary, digits), " "= stars)
    names(printTrms) <- c("P(gamma = 1)", "pi", "dim", " ")
    print(printTrms, digits = digits)
    wrapCat("*:P(gamma = 1)>.25  **:P(gamma = 1)>.5  ***:P(gamma = 1)>.9")
    cat("\n")
  }

  if(printModels) {
    cat(sprintf("\nPosterior model probabilites (inclusion threshold = %g):\n",
      x$thresh))
    modelTable <- {
      #make ("x","")-vectors out of model names
      models <- sapply(names(x$modelTable), function(x) {
        sapply(strsplit(gsub("1", "x", x),""),
          function(y) gsub("0", "", y))
      })

      models <- rbind(round(x$modelTable, 3), models,
        round(cumsum(x$modelTable), 3))
      rownames(models) <- c("prob.:",
        names(x$predvars[!grepl("u(", names(x$predvars), fixed = TRUE)]),
        "cumulative:")
      models <- data.frame(models)

      # show at least 10/at most 20 models, or enough to represent total .9
      # probability mass.
      if(is.null(showModels)) showModels <- min(20, max(10,
        min(which(cumsum(x$modelTable)>=.9))), NCOL(models))

      models <- models[, 1:showModels, drop = F]
      colnames(models) <- 1:NCOL(models)
      models
    }
    print(modelTable)
  }
  cat("\n")
  invisible(x)
}
