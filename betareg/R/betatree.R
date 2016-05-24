###########################
## Beta regression trees ##
###########################

## high-level convenience interface
betatree <- function(formula, partition, data, subset = NULL, na.action = na.omit, 
		    link = "logit", link.phi = "log",
 		    control = betareg.control(), ...)
{
  ## beta regression trees rely on party package
  stopifnot(require("party"))

  ## transform formula
  formula <- Formula(formula)
  stopifnot(length(formula)[1] == 1L & length(formula)[2] >= 1L)
  if(missing(partition)) {
    stopifnot(length(formula)[2] == 3L)
    partition <- formula(formula, lhs = 0, rhs = 3)
  }
  if(length(formula)[2] == 1L) {
    precision <- ~ 1
  } else {
    precision <- formula(formula, lhs = 0, rhs = 2)
    formula <- formula(formula, lhs = 1, rhs = 1)
  }

  ## formula/data/model pre-processing
  betamod <- betaReg(control)
  ff <- modeltools::dpp(betamod, formula,
    other = list(part = partition, precision = precision), 
    data = data, subset = subset, na.action = na.action)

  ## call mob()
  rval <- party::mob(ff, model = betamod,
    link = link, link.phi = link.phi,
    control = party::mob_control(objfun = function(object) - as.vector(logLik(object)), ...))

  ## add class and return
  structure(list(mob = rval), class = "betatree")
}

## convenience plotting
plot.betatree <- function(x, terminal_panel = party::node_bivplot, tnex = 2,
  pval = TRUE, id = TRUE, ...) {
  plot(x$mob, terminal_panel = terminal_panel, tnex = tnex,
    tp_args = list(id = id, ...), ip_args = list(pval = pval, id = id))
}

## hand-crafted "Next()" to bridge to
## un-exported S4 classes "mob"/"BinaryTree", argh!
logLik.betatree <- function(object, ...) logLik(object$mob, ...)
sctest.betatree <- function(x, ...) strucchange::sctest(x$mob, ...)
weights.betatree <- function(object, ...) weights(object$mob, ...)
summary.betatree <- function(object, ...) summary(object$mob, ...)
print.betatree <- function(x, ...) {
  print(x$mob, ...)
  invisible(x)
}

## infrastructure (copied from party)
terminal_nodeIDs <- function(node) {
  if(node$terminal) return(node$nodeID)
  ll <- terminal_nodeIDs(node$left)
  rr <- terminal_nodeIDs(node$right)
  return(c(ll, rr))
}

## parameters for beta-regression trees
coef.betatree <- function (object, node = NULL, ...) 
{
  object <- object$mob
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  rval <- sapply(party::nodes(object, node), function(z) coef(z$model, ...))
  if (!is.null(dim(rval))) {
    rval <- t(rval)
    rownames(rval) <- node
  }
  return(rval)
}



############################################
## S4 StatModel object for betaregression ##
############################################

## StatModel creator function for plug-in to mob()
betaReg <- function(control = betareg.control()) {
  stopifnot(requireNamespace("modeltools"))
  new("StatModel",
    capabilities = new("StatModelCapabilities"),
    name = "beta regression model",
    dpp = modeltools::ModelEnvFormula,
    fit = function(object, weights = NULL, ...) {
      ## extract design and response matrix from the `ModelEnv' object
      y <- as.vector(object@get("responseMatrix"))
      xmat <- object@get("designMatrix")
      zmat <- model.matrix(object@formula$precision, object@get("precision"))
      attr(xmat, "assign") <- attr(zmat, "assign") <- NULL
      names(y) <- rownames(xmat)
    
      ## beta regression fit
      z <- betareg.fit(xmat, y, zmat, weights = weights, control = control, ...)

      ## set class
      class(z) <- c("betaReg", "betareg")

      ## compute empirical estimating functions
      wts <- if(is.null(weights)) 1 else weights
      ## compute y*    
      ystar <- qlogis(y)
      ## compute mu*
      eta <- as.vector(xmat %*% z$coefficients$mean)
      phi_eta <- as.vector(zmat %*% z$coefficients$precision)
      mu <- z$link$mean$linkinv(eta)
      phi <- z$link$precision$linkinv(phi_eta)
      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
      ## compute scores of beta
      ef <- phi * (ystar - mustar) * as.vector(z$link$mean$mu.eta(eta)) * wts * xmat
      ## combine with scores of phi
      if(control$phi) {
        ef <- cbind(ef,
          (mu * (ystar - mustar) + log(1-y) - digamma((1-mu)*phi) + digamma(phi)) *
          as.vector(z$link$precision$mu.eta(phi_eta)) * wts * zmat)
        colnames(ef) <- names(coef(z, phi = control$phi))
      }
      z$estfun <- ef

      ## reconstruct correct formula
      ff <- y ~ .
      ff[[2]] <- object@formula$response[[2]]    
      ff <- as.Formula(update(ff, as.formula(object@formula$input)),
        as.formula(object@formula$precision))
      z$formula <- formula(ff)    
      z$contrasts <- list(mean = attr(xmat, "contrasts"), precision = attr(zmat, "contrasts"))
      z$terms <- list(
        mean = attr(object@get("input"), "terms"),
        precision = attr(object@get("precision"), "terms"),
        full = terms(formula(ff, collapse = TRUE)))
    
      z$addargs <- list(...)
      z$ModelEnv <- object
      z
    }
  )
}

## modify/simply some methods
estfun.betaReg <- function(object, ...) {
    object$estfun
}

print.betaReg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("betaReg fit with coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
    invisible(x)
}

summary.betaReg <- function(object, type = "response", ...)
  summary.betareg(object, type = type, ...)

reweight.betaReg <- function(object, weights, ...) {
     fit <- betaReg(betareg.control(
       phi = object$phi,
       method = object$method,
       maxit = object$control$maxit,
       trace = object$control$trace,
       start = object$coefficients))@fit
     do.call("fit", c(list(object = object$ModelEnv, weights = weights), object$addargs))
}
