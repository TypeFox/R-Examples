# "penfit" object to store the result of a penalized regression
# These objects are meant to be accessed by the user, but not made directly by the user
setClass("penfit", 
  representation(
    penalized = "vector", 
    unpenalized = "vector",
    residuals = "vector",
    fitted = "vector",
    lin.pred = "vector",
    loglik = "numeric",
    penalty = "vector",
    iterations = "numeric",
    converged = "logical",
    model = "character",
    lambda1 = "vector",
    lambda2 = "vector",
    nuisance = "list",
    weights = "vector",
    formula = "list" 
  )
)

# creation method for a penfit object 
.makepenfit <- function(object, unpenalized, model, lambda1, lambda2, fusedl, orthogonalizer, weights, formula) {
  out <- new("penfit")
                                              
  if(!fusedl){
  object$beta <- object$beta / weights}

  beta <- object$beta[unpenalized + seq_len(length(object$beta) - unpenalized)]
  gamma <- object$beta[seq_len(unpenalized)] - drop(orthogonalizer %*% beta)
   
  out@unpenalized <- gamma
  out@penalized <- beta


  out@residuals <- object$fit$residuals
  out@fitted <- object$fit$fitted
  out@lin.pred <- object$fit$lp
                                    
  out@loglik <- if (is.na(object$fit$loglik)) -Inf else object$fit$loglik
  out@penalty <- object$penalty
  
  out@iterations <- object$iter
  
  out@converged <- object$converged
  
  out@model <- model
  
  out@formula <- formula
  
  if ("baseline" %in% names(object$fit$nuisance))
    object$fit$nuisance$baseline <- object$fit$nuisance$baseline()
  out@nuisance <- object$fit$nuisance
  
  out@lambda1 <- lambda1
  out@lambda2 <- lambda2
  out@weights <- weights
  
  out
}

# show method
setMethod("show", "penfit", function(object) {
  cat("Penalized", object@model, "regression object\n")
  if (object@converged) {
    coefs <- unlist(c(object@penalized, object@unpenalized))
    cat(length(coefs), "regression coefficients")
    if (any(coefs == 0)) cat(" of which", sum(coefs != 0), "are non-zero")
    cat("\n\n")
    cat("Loglikelihood =\t", object@loglik, "\n")
    if (any(object@lambda1 > 0))
      if (length(object@lambda1) > 3) 
        cat("L1 penalty =\t", object@penalty[1], "\tat lambda1 = ", object@lambda1[1:3], "...\n")
      else
        cat("L1 penalty =\t", object@penalty[1], "\tat lambda1 = ", object@lambda1, "\n") 
    if (any(object@lambda2 > 0))
      if (length(object@lambda2) > 3) 
        cat("L2 penalty =\t", object@penalty[2], "\tat lambda2 = ", object@lambda2[1:3], "...\n")
      else
        cat("L2 penalty =\t", object@penalty[2], "\tat lambda2 = ", object@lambda2, "\n") 
  } else {
    cat("Model failed to converge\n")
  }
})

# extracts the coefficients
setMethod("coefficients", "penfit", function(object, which = c("nonzero", "all", "penalized", "unpenalized"), standardize = FALSE) {
  which <- match.arg(which)
  nunp <- length(object@unpenalized)
  np <- length(object@penalized)
  whichunp <- switch(which, 
    all =, unpenalized =, nonzero = rep(TRUE,nunp),
    penalized = rep(FALSE, nunp))
  whichp <- switch(which,
    all =, penalized = rep(TRUE, np),
    unpenalized = rep(FALSE, np),
    nonzero = (object@penalized != 0))
  out <- c(object@unpenalized[whichunp], object@penalized[whichp])
  if (standardize) out <- out * object@weights[c(whichunp, whichp)]
  out
})
setMethod("coef", "penfit", function(object, which = c("nonzero", "all", "penalized", "unpenalized"), standardize = FALSE) {
  coefficients(object, which, standardize)
})


# extracts the residuals
setMethod("residuals", "penfit", function(object, ...) {
  object@residuals
})

# extracts the linear predictors
setGeneric("linear.predictors", function(object, ...) standardGeneric("linear.predictors"))
setMethod("linear.predictors", "penfit", function(object, ...) {
  object@lin.pred
})



# extracts the fitted values
setMethod("fitted.values", "penfit", function(object, ...) {
  object@fitted
})
setMethod("fitted", "penfit", function(object, ...) {
  object@fitted
})

# extracts the weights
setMethod("weights", "penfit", function(object, ...) {
  object@weights
})

# extracts the baseline hazard (survival models only)
setGeneric("basesurv", function(fit, centered = TRUE, ...) standardGeneric("basesurv"))
setMethod("basesurv", "penfit", function(fit, centered = TRUE) {
  if (fit@model == "cox") 
    if (centered) {
      meanlp <- mean(linear.predictors(fit))  
      out <- fit@nuisance$baseline
      out@curves <- out@curves^exp(meanlp)
      return(out)
    } else {
      return(fit@nuisance$baseline)
    }
  else
    return(NULL)
})

# calculates a cumulative hazard from a survival curve
setGeneric("basehaz", package = "survival")
setMethod("basehaz", "penfit", function(fit, centered = TRUE) {
  if (fit@model == "cox") {
    bs <- basesurv(fit, centered)
    if (nrow(bs@curves) == 1) rownames(bs@curves) <- "hazard"
    out <- data.frame(-log(t(bs@curves)), time = time(bs), check.names=FALSE)
    return(out)
  } else
    return(NULL)
})

# extracts the penalty
setGeneric("penalty", function(object, ...) standardGeneric("penalty"))
setMethod("penalty", "penfit", function(object, ...) {
  object@penalty
})

# extracts the likelihood
setGeneric("loglik", function(object, ...) standardGeneric("loglik"))
setMethod("loglik", "penfit", function(object, ...) {
  object@loglik
})

# predicts on new data
setMethod("predict", "penfit", function(object, penalized, unpenalized, data) {
                                                        
  # determine defaults
  has.offset <- length(attr(terms(object@formula$unpenalized), 'offset')) != 0
  if (missing(unpenalized)) {
    if (length(object@unpenalized) == 0 && !has.offset)
      unpenalized <- ~0
    else if (!is.null(object@formula$unpenalized))
      unpenalized <- object@formula$unpenalized
    else if (length(object@unpenalized) == 1 && object@model != "cox" && !has.offset)
      unpenalized <- ~1
    else  
      stop("argument \"unpenalized\" is missing.")
  }
  if (missing(data)) data <- NULL

  # coerce unpenalized into a matrix and find the offset and strata terms
  if (is.data.frame(unpenalized) || is.vector(unpenalized)) {
    if (all(sapply(unpenalized, is.numeric))) {
      unpenalized <- as.matrix(unpenalized)
    } else {
      stop("argument \"unpenalized\" could not be coerced into a matrix")
    }
    offset <- 0
  }
  if (is(unpenalized, "formula")) {
    if(is.null(data)) data <- as.data.frame(matrix(,nrow(penalized),0))
    offset <- model.offset(model.frame(unpenalized, data=data))
    unpenalized <- terms(unpenalized, specials='strata')
    if (is.null(offset)) offset <- 0
    # suppress intercept if necessary
    if (object@model == "cox") {
      if (length(attr(unpenalized, "specials")$strata) > 0) {
        strata <- untangle.specials(unpenalized, "strata", 1)
        strata.nrs <- strata$terms                            # indices of the strata variables in the terms object
        strata.nrs2 <- attr(unpenalized, "specials")$strata   # indices of the strata variables in attr(unpenalized, "variables")
        strata <- eval(parse(text=attr(unpenalized, "term.labels")[strata.nrs]),data)
        unpenalized <- unpenalized[-strata.nrs]
      } else strata <- NULL
      attr(unpenalized, "intercept") <- 1
    } 
    unpenalized <- model.matrix(unpenalized, data)
    if (object@model == "cox") 
      unpenalized <- unpenalized[,-1,drop=FALSE]
  }
 
  # coerce penalized into a matrix
  if (missing(penalized)) 
    if (!is.null(object@formula$penalized))
      penalized <- object@formula$penalized
    else
      stop("\"penalized\" argument is missing with no default.")
  if (is.data.frame(penalized) || is.vector(penalized))
    if (all(sapply(penalized, is.numeric))) {
      penalized <- as.matrix(penalized)
    } else {
      stop("argument \"penalized\" could not be coerced into a matrix")
    }
  if (is(penalized, "formula")) {
    has.intercept <- attr(terms(penalized, data=data), "intercept") == 1
    if(is.null(data)) data <- as.data.frame(matrix(,nrow(unpenalized),0))
    oldcontrasts <- unlist(options("contrasts"))
    options(contrasts = c(unordered = "contr.none", ordered = "contr.diff"))
    penalized <- terms(penalized, data=data)
    # suppress intercept
    attr(penalized, "intercept") <- 1
    penalized <- model.matrix(penalized, data)
    if (has.intercept) penalized <- penalized[,-1,drop=FALSE]
    options(contrasts = oldcontrasts)
  }
  
  # find n
  n <- max(nrow(penalized), nrow(unpenalized))
  if (nrow(penalized) != nrow(unpenalized))
    stop("row counts of \"penalized\", \"unpenalized\" and/or \"data\" do not match")

  # check if dimensions and names match
  if (length(object@penalized) != ncol(penalized))
    stop("the dimension of \"penalized\" does not match the fitted model object")
  if (!is.null(names(object@penalized)) && !is.null(colnames(penalized))
    && !all(names(object@penalized) == colnames(penalized)))
    if (setequal(names(object@penalized), colnames(penalized)))
      penalized <- penalized[,names(object@penalized)]
    else
      warning("variable names in \"penalized\" do not match those in the fitted model")
  if (length(object@penalized) != ncol(penalized))
    stop("the dimension of \"unpenalized\" does not match the fitted model object")
  if (!is.null(names(object@unpenalized)) && !is.null(colnames(unpenalized)) &&
    !all(names(object@unpenalized) == colnames(unpenalized)))
    if (setequal(names(object@unpenalized), colnames(unpenalized)))
      unpenalized <- unpenalized[,names(object@unpenalized)]
    else
      warning("variable names in \"unpenalized\" do not match those in the fitted model")

  # find the linear predictors
  lp <- offset + drop(penalized %*% object@penalized) + drop(unpenalized %*% object@unpenalized)
  
  # find the predictions
  predictions <- switch(object@model,
    cox = .coxpredict(lp, object@nuisance, strata),
    linear = .lmpredict(lp, object@nuisance),
    logistic = .logitpredict(lp, object@nuisance),
    poisson = .poissonpredict(lp, object@nuisance)
  )
  
  return(predictions)
})