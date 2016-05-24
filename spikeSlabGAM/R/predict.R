#' @include spikeSlabGAM.R
{}

#' Get the posterior distribution of the linear predictor of a model term
#'
#' @param label (character) one of the terms in \code{model}
#' @param model  a \code{spikeSlabGAM} object
#' @param betaInd (optional) the column indices of the part of the design matrix
#'   for which the linear predictor is to be evaluated.
#' @return a matrix containing the samples of the linear predictor associated
#'   with \code{label}; with attribute \code{'coef'} that contains the posterior
#'   samples of the associated coefficients.
#' @author Fabian Scheipl
#' @export
getPosteriorTerm <- function(label = NULL, model, betaInd = NULL) {
  stopifnot(any(c(!is.null(label), !is.null(betaInd))),
    class(model)=="spikeSlabGAM")
  if(is.null(betaInd)) {
    betaInd <- which(model$model$groupIndicatorsOrig == label)
  }

  beta <- t(do.call(rbind,
    lapply(model$samples$beta, "[",
      i = 1:model$mcmc$chainLength, j = betaInd, drop = F)))
  return(structure(model$X[, betaInd, drop = F] %*% beta, "coef"= beta))
}

#' Get summaries of the posterior (predictive) distribution of the linear
#' predictor of a model term
#'
#' @param label (character) the label of one of the terms in \code{model}.
#' @param model  a \code{spikeSlabGAM} object
#' @param newdata \code{data.frame} on which to evaluate \code{label}. Defaults
#'   to NULL, in which case the term is evaluated on the original data.
#' @param aggregate (function) a summary statistic that is applied over the
#'   mcmc-samples of the linear predictor. Defaults to \code{mean}.
#' @param quantiles (numeric) a vector of quantiles for borders of credible
#'   regions of the linear predictor. Defaults to 10 and 90 percent quantiles,
#'   i.e. a (point-wise) 80 percent credible region.
#' @param returnData  should the relevant original variables be included in the
#'   returned \code{data.frame}? Defaults to TRUE.
#' @return A \code{data.frame} that contains the relevant variables from newdata
#'   (if \code{returnData} is TRUE), the \code{aggregate}-summary and the
#'   requested \code{quantiles}.
#' @author Fabian Scheipl
#' @export
evalTerm <- function(label, model, newdata = NULL, aggregate = mean,
  quantiles = c(.1,.9), returnData = TRUE) {

  stopifnot(is.null(aggregate)|is.function(aggregate),
    is.null(quantiles)|all(quantiles >= 0 & quantiles <= 1),
    class(model)=="spikeSlabGAM")

  # find main effect(s) and beta
  trms <- strsplit(label, ":")[[1]]
  vars <- sapply(trms, function(x) gsub("\\)", "",
    gsub("([[:graph:]]+\\()","", x)))
  betaInd <- which(names(model$model$groupIndicatorsOrig)== label)

  # get posterior samples of term
  fhat <-	getPosteriorTerm(label, model, betaInd)
  if(!is.null(newdata)) {
    # make (uncentered, high-rank) design(s) for newdata
    newX <- vector(length(trms), mode ="list")
    names(newX) <- trms
    for(trm in trms) {
      newX[[trm]] <- do.call(model$predvars[[trm]]$makeNewX,
        list(xnew = newdata[, vars[trm]]))
    }
    if(length(newX)>1) {
      #make interaction design
      X <- matrix(1, nrow = nrow(newX[[1]]))
      for(i in 1:length(newX)) {
        X <- matrix(apply(X, 2, function(x) {
          x * newX[[i]]}),
          nrow = nrow(X), ncol = NCOL(X)* NCOL(newX[[i]]))
      }
    } else {
      X <- newX[[1]]
    }
    coef <- if(!is.null(model$predvars[[label]]$qrB)) {
      # posterior of coefficients for un-orthogonalized interactions: idea: 1. B
      # is the design for the interaction BEFORE low-rank approx. (saved in
      # predvars) 2. fhat is the posterior distribution of B %*% coef (sampled via
      # B_loRank %*% coef_loRank) --> get coefs for original design matrix B via
      # qr(B) : R beta = Q' fhat
      tmp <- try(qr.coef(model$predvars[[label]]$qrB, fhat), silent = TRUE)
      if(class(tmp)=="try-error") {
        tmp <- (ginv(qr.R(model$predvars[[label]]$qrB)) %*%
            t(qr.Q(model$predvars[[label]]$qrB)))%*% fhat
        ## FIXME: use backsolve? or do we need ginv for stability
      }
      tmp
    } else {
      attr(fhat, "coef")
    }
    if(any(is.na(coef))) {
      ##FIXME
      warning("Recalculating original coefficients for ", label,
        " seems to be problematic. Predictions may be unstable." )
      rm <- is.na(coef[, 1])
      #coef <- coef[!rm, , drop = F]
      #X <- X[, !rm, drop = F]
      coef[rm, ] <- 0
    }
    #posterior of term evaluated on newdata:
    fhat <- X %*% coef
    data <- newdata[, vars, drop = F]
  } else {
    data <- model$data[, vars, drop = F]
  }
  pred <- summarizeF(fhat, aggregate, quantiles)
  if(returnData) {
    return(cbind(data, pred))
  } else {
    return(pred)
  }
}

#' Obtain posterior predictive/credible intervals from a spike-and-slab model
#'
#' @param object a \code{spikeSlabGAM} model
#' @param newdata an optional \code{data.frame} on which to evaluate
#'   predictions. Defaults to NULL, in which case the fitted values for the
#'   original data are returned
#' @param type 	the type of prediction required. The default is on the scale of
#'   response, the alternative \code{'link'} is on the scale of the linear
#'   predictor. Specifying \code{'terms'} returns a list giving the linear
#'   predictors of the terms specified in \code{terms}.
#' @param terms an optional character vector of term labels or variable names
#'   for which to return fits/predictions/credible regions. If variable names
#'   instead of term labels are supplied, the function returns
#'   predictions/estimates for \emph{all} terms associated with these variables,
#'   i.e. their main effects (usually both linear and smooth for numeric
#'   covariates) and all interaction terms they are involved in.
#' @param aggregate (function) the summary statistic of the posterior predictive
#'   of the linear predictor. Defaults to \code{mean}.
#' @param quantiles (numeric) an optional vector of quantiles for borders of
#'   credible regions of the returned values. Defaults to NULL.
#' @param addIntercept include global intercept term in prediction/estimate?
#'   Defaults to TRUE if \code{terms = NULL}.
#' @param ... 	arguments passed from or to other methods (not used)
#' @return  If \code{type ="terms"}, a list of \code{data.frame}s containing the
#'   requested pointwise summary statistics for the supplied terms (use e.g.
#'   \code{\link{Reduce}("+", ...)} to get row-wise sums of the list-entries).
#'   Otherwise, a \code{data.frame} containing the requested pointwise summary
#'   statistics of the posterior predictive of the linear predictor
#'   (\code{type ="link"}) or the conditional expectation of the response
#'   (\code{type ="response"}) is returned.
#' @export
#' @author Fabian Scheipl
predict.spikeSlabGAM <- function(object, newdata = NULL,
  type = c("response", "link", "terms"), terms = NULL,
  aggregate = mean, quantiles = NULL, addIntercept = is.null(terms), ...) {

  addIntercept <- eval(addIntercept)
  type <- match.arg(type)

  if(!is.null(terms)) type <- "terms"
  linkinv <- switch(as.character(object$family),
    "0"= I,
    "1"= plogis,
    "2"= exp)

  if(is.null(newdata)) {
    if(type!="terms") {
      fhat <- getPosteriorTerm(model = object, betaInd = 1:object$model$q)
      res <- summarizeF(fhat, aggregate, quantiles)
      if(type =="response") {
        res <- linkinv(res)
      }
      return(res)
    } else {
      terms <- if(is.null(terms)) {
        names(object$predvars)
      } else {
        #check for special terms
        specials <- sapply(terms,
          function(x) any(mapply(grepl, x = x,
            pattern = environment(object$formula)$specials)))
        #find out terms that involve variables mentioned in non-special "terms"
        which <- unique(unlist(sapply(terms[!specials], grep,
          x = names(object$predvars), fixed = TRUE)))
        unique(c(terms[specials], names(object$predvars)[sort(which)]))
      }
      res <- lapply(terms, function(x) {
        tmp <- getPosteriorTerm(label = x, model = object)
        summarizeF(tmp, aggregate = aggregate, quantiles = quantiles)
      })
      names(res) <- terms
      if(attr(terms(object$formula),"intercept")!= 0 & addIntercept) {
        int <- getPosteriorTerm(model = object, betaInd = 1)
        int <- summarizeF(int, aggregate = aggregate, quantiles = quantiles)
        res[["(Intercept)"]] <- int
      }
      return(res)
    }
  } else {
    terms <- if(type != "terms" | is.null(terms)) {
      names(object$predvars)
    } else {
      #check for special terms
      specials <- sapply(terms,
        function(x) any(mapply(grepl, x = x,
          pattern = environment(object$formula)$specials)))
      #find out terms that involve variables mentioned in non-special "terms"
      which <- unique(unlist(sapply(terms, grep,
        x = names(object$predvars), fixed = TRUE)))
      names(object$predvars)[sort(which)]
    }
    res <- lapply(terms, evalTerm, model = object, newdata = newdata,
      aggregate = aggregate, quantiles = quantiles, returnData = FALSE)
    if(type =="terms") {
      names(res) <- terms
      if(attr(terms(object$formula),"intercept")!= 0 & addIntercept) {
        int <- getPosteriorTerm(model = object, betaInd = 1)
        int <- summarizeF(int, aggregate = aggregate, quantiles = quantiles)
        res[["(Intercept)"]] <- suppressWarnings(matrix(int, nrow = nrow(newdata),
          ncol = ncol(int), byrow = T))
      }
      return(res)
    } else {
      if(attr(terms(object$formula),"intercept")!= 0 & addIntercept) {
        int <- getPosteriorTerm(model = object, betaInd = 1)
        int <- summarizeF(int, aggregate = aggregate, quantiles = quantiles)
        res[["(Intercept)"]] <- suppressWarnings(matrix(int, nrow = nrow(newdata),
          ncol = ncol(int), byrow = T))
      }
      res <- Reduce("+", res)
      if(type =="response") {
        res <- linkinv(res)
      }
      return(res)
    }
  }
}
