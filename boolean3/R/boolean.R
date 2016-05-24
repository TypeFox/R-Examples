### ----------------------------------------------------------------------------
### This file is part of boolean3
###
### Copyright (C) 2011--2014 Jason W. Morgan <morgan.746@osu.edu>
###
### boolean3 represents a substantial re-write of the original boolean package
### developed by Bear Braumoeller, Ben Goodrich, and Jacob Kline. This version
### was developed under the direction of Bear Braumoeller and with support from
### The Ohio State University's College of Social and Behavioral Sciences.
###
### boolean3 is free software: you can redistribute it and/or modify it under
### the terms of the GNU General Public License as published by the Free
### Software Foundation, either version 3 of the License, or (at your option)
### any later version.
###
### This program is distributed in the hope that it will be useful, but WITHOUT
### ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
### FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
### more details.
###
### You should have received a copy of the GNU General Public License along with
### this program.  If not, see <http://www.gnu.org/licenses/>.
###
### ----------------------------------------------------------------------------

##' Construct a boolean object that can then be fit with \code{boolean}.
##'
##' \code{boolprep} sets up a boolean model object that can then be fit with
##' \code{\link{boolean}}. A properly specified model (contained in \code{...})
##' will contain at least three components. The first component must specify the
##' boolean logic to be employed. For instance, the \code{y ~ (a | b)} formula
##' would indicate a logical \code{or} between the \code{a} and \code{b}
##' submodels, while \code{y ~ (a & b)} would indicate a logical
##' \code{and}. \code{y} is the name of the response variable of
##' interest. Logical operators can be nested; e.g., \code{y ~ (a | (b & c))} is
##' valid. The second and third components are submodels and are specified as
##' usual: \code{a ~ x1 + x2} and \code{b ~ x3 + x4 + x5}, where are the
##' \code{x}-variables are covariates. \code{a}, \code{b}, and \code{c} are
##' labels indicating the submodel position in the boolean specification.
##' @title Specify Boolean Model
##' @param ... formula specification for boolean model. See the details and
##' examples sections below.
##' @param data \code{data.frame} containing the data to be used in the model.
##' @param subset select subset of data. See \code{\link{lm}} for details.
##' @param weights specify model weights (not implemented).
##' @param na.action set \code{na.action}. See \code{\link{lm}} for details.
##' @param offset specify an offset. Offsets are not implemented and this
##' parameter is simply ignored.
##' @param family a model \code{family} to use for estimation. This can be
##' comprised of a list of link functions when the desire is to have model
##' components with different links. \code{binomial} is the only family
##' supported at this time. See \code{\link{family}} for more details.
##' @return \code{boolprep} returns a \code{\link{boolean-class}} object
##' containing the model components needed to estimate a boolean model.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @references Braumoeller, Bear F. (2003) ``Causal Complexity and the Study
##' of Politics.'' \emph{Political Analysis} 11(3): 209--233.
##' @seealso See \code{\link{boolean}} for model estimation.
##' @examples
##' ## Generate some fake data
##' require(mvtnorm)
##' set.seed(12345)
##' N  <- 2000
##' Df <- cbind(1, rmvnorm(N, mean=rep(0, 5)))
##'
##' ## Set coefficients
##' beta.a <- c(-2.00, 0.33, 0.66, 1.00)
##' beta.b <- c(0.00, 1.50, -0.25)
##'
##' ## Generate path probabilities following a normal model.
##' y.a <- as.vector(pnorm(tcrossprod(beta.a, Df[, 1:4])))
##' y.b <- as.vector(pnorm(tcrossprod(beta.b, Df[, c(1, 5, 6)])))
##'
##' ## AND and OR-model
##' or <- function(x, y) { x + y - x * y }
##' and <- function(x, y) { x * y }
##' y.star.OR  <- or(y.a, y.b)
##' y.star.AND <- and(y.a, y.b)
##'
##' ## Observed responses
##' y.OR <- rbinom(N, 1, y.star.OR)
##' y.AND <- rbinom(N, 1, y.star.AND)
##'
##' ## Set up data.frame for estimation
##' Df <- cbind(1, Df)
##' Df <- as.data.frame(Df)
##' Df[,1] <- y.OR
##' Df[,2] <- y.AND
##' names(Df) <- c("y.OR", "y.AND", "x1", "x2", "x3", "x4", "x5")
##' 
##' ## Before estimating, boolean models need to be specified using the
##' ## boolprep function.
##' 
##' ## OR model, specified to use a probit link for each causal path. This
##' ## matches the data generating process above.
##' mod.OR <- boolprep(y.OR ~ (a | b), a ~ x1 + x2 + x3, b ~ x4 + x5,
##'                    data = Df, family=list(binomial("probit")))
##' @export
boolprep <- function(..., data, subset, weights, na.action, offset,
                     family=list(binomial(link="logit"))) {

    call <- match.call(expand.dots=TRUE) # Get model call
    frame <- build_frame(call)           # Build the model frame
    model <- build_model(call)           # Build model specification
    links <- make.links(family, model)   # Make link function(s)

    ## Set up object and return.
    obj <-
        structure(list
                  (call = call,
                   link = links,
                   model = model,
                   N = nrow(frame$frame[[1]]),         # Observations
                   k = sum(sapply(frame$frame, ncol)), # Number of coefficients
                   coef.labels = set_labels(frame$frame), # Coefficient labels
                   coef.idx = set_coef(frame$frame),      # Coefficient index
                   response = frame$response,             # Extract response
                   frame = frame$frame),                  # Model frame
                  class="boolean")

    ## Set up scobit lambda (shape) parameter(s) if necessary. coef.idx and k
    ## need updated if this is the case. CURRENTLY NOT AVAILABLE.
    if ("scobit" %in% mapply(function(x) x$name, links)) {
        ## obj <- scobit_setup(obj)
        stop("scobit link function is not supported.")
    }

    ## Build calc_lik closure.
    obj$calc_lik <- make_lik(obj)
  
    obj
}

##' Performs a fit of a boolean model as specified by \code{\link{boolprep}}.
##'
##' \code{boolean} performs a fit of a boolean model as specified by
##' \code{\link{boolprep}}.
##' @title Fit a Boolean Model
##' @param obj boolean model object as produced by \code{\link{boolprep}}.
##' @param method string (or string vector) specifying the method(s) of
##' estimation. The specified method(s) should be one of those available from
##' the \code{\link{optimx}} or \code{\link{optim}} functions. A genetic
##' algorithm is available from \code{\link{genoud}}
##' (\code{method="genoud"}). \code{method} defaults to \code{"nlminb"}.
##' @param start numeric vector specifying starting values for each parameter in
##' the model. Must have a length equal to the number of parameters being
##' estimated. Defaults to \code{NULL}, which instructs \code{boolean} to
##' estimate ``sensible'' starting values (currently the coefficient values
##' estimated from a \code{\link{glm}} model).
##' @param ... additional arguments to be passed on to optimizers. Each
##' optimizer provides numerous optional parameters to help improve estimation
##' results. See the provided examples and the documentation for the estimation
##' method of interest.
##' @return A boolean model object containing the fit results (detailed results
##' available in the \code{model.fit} slot).
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @references Braumoeller, Bear F. (2003) ``Causal Complexity and the Study of
##' Politics.'' \emph{Political Analysis} 11(3): 209--233.
##' @seealso See \code{\link{boolprep}} for model setup, \code{\link{boolean}}
##' the \code{snow} package for details on clustering (useful when using
##' \code{\link{genoud}}). See \code{\link{optimx}}, \code{\link{optim}}, and
##' \code{\link{genoud}} for detailed documentation on the estimation methods
##' available.
##' @examples
##' ## Generate some fake data
##' require(mvtnorm)
##' set.seed(12345)
##' N  <- 2000
##' Df <- cbind(1, rmvnorm(N, mean=rep(0, 5)))
##'
##' ## Set coefficients
##' beta.a <- c(-2.00, 0.33, 0.66, 1.00)
##' beta.b <- c(0.00, 1.50, -0.25)
##'
##' ## Generate path probabilities following a normal model.
##' y.a <- as.vector(pnorm(tcrossprod(beta.a, Df[, 1:4])))
##' y.b <- as.vector(pnorm(tcrossprod(beta.b, Df[, c(1, 5, 6)])))
##'
##' ## AND and OR-model
##' or <- function(x, y) { x + y - x * y }
##' and <- function(x, y) { x * y }
##' y.star.OR  <- or(y.a, y.b)
##' y.star.AND <- and(y.a, y.b)
##'
##' ## Observed responses
##' y.OR <- rbinom(N, 1, y.star.OR)
##' y.AND <- rbinom(N, 1, y.star.AND)
##'
##' ## Set up data.frame for estimation
##' Df <- cbind(1, Df)
##' Df <- as.data.frame(Df)
##' Df[,1] <- y.OR
##' Df[,2] <- y.AND
##' names(Df) <- c("y.OR", "y.AND", "x1", "x2", "x3", "x4", "x5")
##' 
##' ## Before estimating, boolean models need to be specified using the
##' ## boolprep function.
##' 
##' ## OR model, specified to use a probit link for each causal path. This
##' ## matches the data generating process above.
##' mod.OR <- boolprep(y.OR ~ (a | b), a ~ x1 + x2 + x3, b ~ x4 + x5,
##'                    data = Df, family=list(binomial("probit")))
##' 
##' ## Fit a model using the nlminb optimizer (the default). Verbose output is
##' ## requested by specifying trace=1 as a control parameter.
##' (fit <- boolean(mod.OR, method="nlminb", control=list(trace=1)))
##' 
##' ## Multiple optimizers can be specified in a single call to boolean. Here 
##' ## we fit with the nlm and nlminb optimizers.
##' (fit1 <- boolean(mod.OR, method=c("nlm", "nlminb")))
##'
##' ## The summary function will report the detailed results for each
##' ## optimization method. 
##' summary(fit1)
##' @export
boolean <- function(obj, method="nlminb", start=NULL, ...) {

    ## Verify that obj is a boolean object.
    if (!inherits(obj, "boolean"))
        stop("boolean requires a boolean object as outputed by boolprep", .call=FALSE)

    ## Assign the estimation method if it doesn't exist.
    if(is.null(obj$method))
        obj$method <- method

    ## Set the starting values if the are not provided.
    if (is.null(start))
        start <- set_start(obj)
    else if (length(start) != obj$k)
        stop("wrong number of starting values specified", .call=FALSE)

    ## Fit the model given the selected optimization method.
    if (length(obj$method) == 1) {
        fit <- switch(as.character(obj$method),
                      "genoud"= fit_genoud(obj, start=start, ...),
                      ## "mcmc"  = fit_mcmc(obj, start=start, ...),
                      "SANN"  = fit_optim(obj, start=start, ...),
                      fit_optimx(obj, start=start, ...))
    }
    else {
        fit <- fit_optimx(obj, start=start, ...)
    }
  
    ## Set the class of the fit object, add it to the boolean object, and return.
    class(fit) <- c("boolfit", class(fit))
    obj$model.fit <- fit
    obj
}


### ----------------------------------------------------------------------------
### Utilities
### ----------------------------------------------------------------------------

make.links <- function(family, model) {
  ## Make link function(s). Family must be a list of length 1 or equal to the
  ## number of submodels. If the latter is the case, the family list must have
  ## the same item names as in model$submods.

  ## All in family have to be of the "family" class.
  if (!(all(Map(class, family) == "family")))
    stop("all elements of 'family' must be of a valid 'family' class", .call=FALSE)

  ## Family has to be of length 1 or the same length as model$submods.
  if (! (length(family) %in% c(1, length(model$submods))))
    stop("incorrect number of model families specified in 'family'", .call=FALSE)
  
  ## If more than one family is specified, their names need to match those in
  ## submods.
  if (length(family) == length(model$submods)) {
    if (!all(names(family) %in% names(model$submods)))
        stop("submodel identifiers in 'family' do not match those of the specified model", .call=FALSE)
  }

  ## If one family specified, replicate it to match the length of the submods.
  if (length(family) == 1) {
    family <- rep(list(family[[1]]), length(model$submods))
    names(family) <- names(model$submods)
  }

  ## Iterate through family list, creating requested link functions.
  Map(function(x) {bool.link(x$link)}, family)  
}

get_rhs <- function(form) {
  ## Get the right hand side of a formula object.
  as.formula(paste("~", as.character(form[3])))
}

get_lhs <- function(form) {
  ## Get the response variable from the formula.
  as.character(form[2])
}

set_labels <- function(frame) {
  ## Set coefficient label names.
  labels <- c()
  for (m in names(frame)) {
    lbl <- paste(colnames(frame[[m]]), m, sep = "_")
    labels <- c(labels, lbl)
  }
  labels
}

build_modframe <- function(formulas) {
  ## Build a single model specification based on the formula inputs.

  forms <- as.character(as.list(formulas)[-1])
  
  ## Check that the first formula is a boolean-type formula. stop() if not.
  if (!isTRUE(grep("[&|]", forms[[1]]) == 1))
      stop("first formula does not appear to be a boolean specification", .call=FALSE)

  ## Check that there are at least three formula entries. stop() if not.
  if (length(forms) < 3)
      stop("insufficient number of formulas", .call=FALSE)
  
  ## Replace boolean operators with `+'. All parentheses to "". Then convert
  ## all formula strings to formula objects.
  forms <- as.list(forms)
  forms[[1]] <- gsub("[&|]", "+", gsub("[()]", "", forms[[1]]))
  for (i in 1:length(forms)) {
    f <- forms[[i]]
    forms[[i]] <- as.formula(forms[[i]])
  }

  ## Verify that all formulas have a response == 1; i.e., they have to be
  ## two-sided formulas. stop() if not.
  for (i in 1:length(forms)) {
    r <- attr(terms(forms[[i]]), "response")
    if (r != 1)
        stop("all formulas must specify a single response variable", .call=FALSE)
  }
  
  ## Check that all components specified in the boolean model exist as
  ## responses in the other submodels. stop() if not.
  v <- rownames(attr(terms(forms[[1]]), "factors"))[-1]
  for (i in 2:length(forms)) {
    r <- rownames(attr(terms(forms[[i]]), "factors"))[1]
    if (!(r %in% v))
        stop("component of the boolean formula not specified", .call=FALSE)
  }

  ## Merge the models into a single simple specification. Variables in the
  ## boolean model are *not* retained (they are latent and so don't actually
  ## exist in the data).
  y <- rownames(attr(terms(forms[[1]]), "factors"))[1]
  model <- c()
  for (i in 2:length(forms)) {
    r <- rownames(attr(terms(forms[[i]]), "factors"))[-1]
    model <- c(model, r)
  }

  model <- paste(y, "~", paste(model, collapse = " + "))
  formula(model)
}

build_frame <- function(mod.call) {
  ## Builds a model.frame from all components included in the call.

  formulas <- mod.call[names(mod.call) %in% ""]
  formula  <- build_modframe(formulas)
  
  ## Get the components needed from the model call.
  m <- match(c("data", "subset", "weights", "na.action", "offset"),
             names(mod.call), 0)

  ## Build the overall model frame. (This bit of code was adapted from lme4.)
  ## This model frame will then be used to construct the submodel frames. This
  ## first full matrix is constructed to assure that there aren't any
  ## mis-matches in obs due to missing data, subsetting, etc.
  frame <- list()
  full.matrix <- mod.call[c(1, m)]           # get the data.frame name.
  full.matrix$formula <- formula             # add the full formula obj.
  full.matrix$drop.unused.levels <- TRUE     # says to drop factor levels
  full.matrix[[1]] <- as.name("model.frame") # renames function to model.frame
  full.matrix <- eval(full.matrix)

  ## Loop through the sub-formulas, constructing model frames as specified by
  ## the user. This guarantees that intercepts can be dropped, for example.
  for (m in 3:length(formulas)) {
    mname <- get_lhs(formulas[[m]])
    rhs <- get_rhs(formulas[[m]])
    frame[[mname]] <- model.matrix(rhs, data = full.matrix)
  }

  list(response = model.response(full.matrix), frame = frame)
}

build_model <- function(mod.call) {
  ## Builds the model, in boolean form, to be estimated.

  ## Set up the container for the models.
  model <- list(submods=list())
  
  ## Get the formulas.
  formulas <- mod.call[names(mod.call) %in% ""]
  forms <- as.list(as.character(as.list(formulas)[-1]))

  ## Set the boolean model, adding the boolean operators and parsing it to
  ## prepare for evaluation.
  model$boolmod <- gsub("[&]", "%&%", gsub("[|]", "%|%", forms[[1]]))
  model$boolmod <- parse(text = as.character(formula(model$boolmod)[3]))

  ## Set the submodels.
  for (i in 2:length(forms)) {
    f <- as.formula(forms[[i]])
    n <- rownames(attr(terms(f), "factors"))[1]
    model$submods[[n]] <- attr(terms(f), "term.labels")
  }
  
  model
}

set_coef <- function(frame) {
  ## Returns a list of index vectors for the coefficient vector for each of
  ## the submodels specified in the frame object.
  coef.idx <- list()

  ## Get the total number of coefficients.
  b <- rep(0, sum(sapply(frame, ncol)))
  i <- 1
  
  for (f in names(frame)) {
    k <- ncol(frame[[f]])
    ##! k <- nrow(frame[[f]])
    c <- b
    c[i:(i+k-1)] <- 1
    coef.idx[[f]] <- as.logical(c)
    i <- i+k
  }
  coef.idx
}

scobit_setup <- function(obj) {
  ## Set up scobit lambda (shape) parameter(s) if necessary. coef.idx,
  ## coef.labels, and k need updated. Returns a modified obj.
  l <- (mapply(function(x) x$name, obj$link) == "scobit") + 0
  names(l) <- paste("lambda", names(l), sep="_")
  obj$scobit.lambda <- l[l==1]
  obj$k <- obj$k + sum(obj$scobit.lambda) # Update k
  obj$coef.idx <- Map(function(x) c(x, rep(FALSE, length(obj$scobit.lambda))),
                      obj$coef.idx)
  obj$coef.labels <- append(obj$coef.labels, names(obj$scobit.lambda))
  obj
}
