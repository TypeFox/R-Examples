#' Multinomial Logistic Regression for Dependent Variables with Unordered Categorical Values
#'
#' Vignette: \url{http://docs.zeligproject.org/en/latest/zeligchoice-mlogit.html}
#' @import methods
#' @export Zelig-bprobit
#' @exportClass Zelig-bprobit

zmlogit <- setRefClass("Zelig-mlogit",
                          contains = "Zelig",
                          field = list(family = "ANY",
                                       linkinv = "function"
                          ))

zmlogit$methods(
  initialize = function() {
    callSuper()
    .self$name <- "mlogit"
    .self$description <- "Multinomial Logistic Regression for Dependent Variables with Unordered Categorical Values"
    .self$fn <- quote(VGAM::vglm)
    .self$authors <- "Matthew Owen, Olivia Lau, Kosuke Imai, Gary King"
    .self$packageauthors <- "Thomas W. Yee"
    .self$year <- 2007
    .self$category <- "multinomial"
    .self$family <- "multinomial"
    .self$wrapper <- "mlogit"
    .self$vignette.url <- "http://docs.zeligproject.org/en/latest/zeligchoice-mlogit.html"
  }
)

zmlogit$methods(
  zelig = function(formula, data, ..., weights = NULL, by = NULL) {
    .self$zelig.call <- match.call(expand.dots = TRUE)
    .self$model.call <- match.call(expand.dots = TRUE)
    .self$model.call$family <- .self$family
    callSuper(formula = formula, data = data, ..., weights = NULL, by = by)
  }
)

zmlogit$methods(
  param = function(z.out) {
    return(mvrnorm(.self$num, coef(z.out), vcov(z.out)))
  }
)

zmlogit$methods(
  # From ZeligChoice 4
  qi = function(simparam, mm) {
    fitted <- .self$zelig.out$z.out[[1]]
    # get constraints from fitted model
    constraints <- fitted@constraints
    coef <- simparam
    ndim <- ncol(fitted@y) - 1
    all.coef <- NULL
    v <- construct.v(constraints, ndim)
    # put all indexed lists in the appropriate section
    for (i in 1:ndim)
      all.coef <- c(all.coef, list(coef[, v[[i]]]))
#     cnames <- ynames <-  if (is.null(colnames(fitted@y))) {1:(ndim + 1)} else colnames(fitted@y)
    if (is.null(colnames(fitted@y))) {
      cnames <- 1:(ndim + 1)
    } else
        cnames <- colnames(fitted@y)
    ynames <- cnames
    cnames <- paste("Pr(Y=", cnames, ")", sep = "")
    ev <- ev.mlogit(fitted, constraints, all.coef, mm, ndim, cnames)
    pv <- pv.mlogit(fitted, ev) #, ynames)
    return(list(ev = ev, pv = pv))
  }
)


#' Split Names of Vectors into N-vectors
#' This function is used to organize how variables are spread
#' across the list of formulas
#' @usage construct.v(constraints, ndim)
#' @param constraints a constraints object
#' @param ndim an integer specifying the number of dimensions
#' @return a list of character-vectors
construct.v <- function(constraints, ndim) {
  v <- rep(list(NULL), ndim)
  names <- names(constraints)
  for (i in 1:length(constraints)) {
    cm <- constraints[[i]]
    for (j in 1:ndim) {
      if (sum(cm[j, ]) == 1) {
        v[[j]] <- if (ncol(cm) == 1)
          c(v[[j]], names[i])
        else
          c(v[[j]], paste(names[i], ':', j, sep=""))
      }
    }
  }
  return(v)
}


#' Simulate Expected Value for Multinomial Logit
#' @usage ev.mlogit(fitted, constraints, all.coef, x, ndim, cnames)
#' @param fitted a fitted model object
#' @param constraints a constraints object
#' @param all.coef all the coeficients
#' @param x a setx object
#' @param ndim an integer specifying the number of dimensions
#' @param cnames a character-vector specifying the names of the columns
#' @return a matrix of simulated values
ev.mlogit <- function (fitted, constraints, all.coef, x, ndim, cnames) {
  if (is.null(x))
    return(NA)
  linkinv <- fitted@family@linkinv
  xm <- rep(list(NULL), ndim)
  sim.eta <- NULL
  x <- as.matrix(x)
  for (i in 1:length(constraints)) {
    for (j in 1:ndim)
      if (sum(constraints[[i]][j,] ) == 1)
        xm[[j]] <- c(xm[[j]], x[, names(constraints)[i]])
  }
  for (i in 1:ndim)
    sim.eta <- cbind(sim.eta, all.coef[[i]] %*% as.matrix(xm[[i]]))
  ev <- linkinv(sim.eta)
  colnames(ev) <- cnames
  return(ev)
}

#' Simulate Predicted Values
#' @usage pv.mlogit(fitted, ev)
#' @param fitted a fitted model object
#' @param ev the simulated expected values
#' @return a vector of simulated values
pv.mlogit <- function (fitted, ev){ #, ynames) {
  if (all(is.na(ev)))
    return(NA)
  # initialize predicted values and a matrix
  pv <- NULL
  Ipv <- sim.cut <- matrix(NA, nrow = nrow(ev), ncol(ev))
  k <- ncol(ev)
  colnames(Ipv) <- colnames(sim.cut) <- colnames(ev)
  sim.cut[, 1] <- ev[, 1]
  for (j in 2:k)
    sim.cut[, j] <- sim.cut[ , j - 1] + ev[, j]
  tmp <- runif(nrow(ev), min = 0, max = 1)
  for (j in 1:k)
    Ipv[, j] <- tmp > sim.cut[, j]
  for (j in 1:nrow(Ipv))
    pv[j] <- 1 + sum(Ipv[j, ])
  pv <- factor(pv, ordered = FALSE)
  pv.matrix <- matrix(pv, nrow = dim(ev)[1])
  levels(pv.matrix) <- levels(pv)
  return(pv.matrix)
}

zmlogit$methods(
  mcfun = function(x, b0=-0.5, b1=0.5, b2=-1, b3=1, ..., sim=TRUE){
    mu1 <- b0 + b1 * x
    mu2 <- b2 + b3 * x

    if(sim){
      n.sim = length(x)
      y.star.1 <- exp( rlogis(n = n.sim, location = mu1, scale = 1) ) # latent continuous y
      y.star.2 <- exp( rlogis(n = n.sim, location = mu2, scale = 1) ) # latent continuous y
      pi1 <- y.star.1/(1 + y.star.1 + y.star.2)
      pi2 <- y.star.2/(1 + y.star.1 + y.star.2)
      pi3 <- 1 - pi1 - pi2

      y.draw <- runif(n=n.sim)
      y.obs <- 1 + as.numeric(y.draw>pi1) + as.numeric(y.draw>(pi1 + pi2))
      return(as.factor(y.obs))
    }else{
      pi1.hat <- exp(mu1)/(1 + exp(mu1) + exp(mu2))
      pi2.hat <- exp(mu2)/(1 + exp(mu1) + exp(mu2))
      pi3.hat <- 1 - pi1.hat - pi2.hat
      
      y.obs.hat <- pi1.hat*1 + pi2.hat*2 + pi3.hat*3    # This is the expectation the MC test will check, although it is not substantively meaningful for factor dep. var.
      return(y.obs.hat)
    }
  }
)
