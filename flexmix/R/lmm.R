setClass("FLXcomponentlmm",
         representation(random="list"),
         contains = "FLXcomponent")

setClass("FLXMRlmm",
         representation(family = "character",
                        random = "formula",
                        group = "factor",
                        z = "matrix",
                        which = "ANY"),
         contains = "FLXMR")

setClass("FLXMRlmmfix",
         contains = "FLXMRlmm")

FLXMRlmm <- function(formula = . ~ ., random, lm.fit = c("lm.wfit", "smooth.spline"),
                     varFix = c(Random = FALSE, Residual = FALSE), ...)
{
  family <- "gaussian"
  lm.fit <- match.arg(lm.fit)
  if (length(varFix) != 2 || is.null(names(varFix)) || any(is.na(pmatch(names(varFix), c("Random", "Residual"))))) 
    stop("varFix has to be a named vector of length two")
  else names(varFix) <- c("Random", "Residual")[pmatch(names(varFix), c("Random", "Residual"))]
  random <- if (length(random) == 3) random else formula(paste(".", paste(deparse(random), collapse = "")))
  object <- new("FLXMRlmm", formula = formula, random = random, 
                weighted = TRUE, family = family, name = "FLXMRlmm:gaussian")
  if (any(varFix)) object <- new("FLXMRlmmfix", object)
  object@preproc.y <- function(x){
    if (ncol(x) > 1)
      stop(paste("y must be univariate"))
    x
  }
  if (lm.fit == "smooth.spline") {
    object@preproc.x <- function(x){
      if (ncol(x) > 1)
        stop(paste("x must be univariate"))
      x
    }
  } 
  add <- function(x) Reduce("+", x)
  lmm.wfit <- function(x, y, w, z, which, random) {
    effect <- lapply(seq_along(which), function(i) z[[which[i]]] %*% random$beta[[i]])
    W <- rep(w, sapply(x, nrow))
    X <- do.call("rbind", x)
    Y <- do.call("rbind", y)
    Effect <- do.call("rbind", effect)
    fit <- get(lm.fit)(X, Y - Effect, W, ...)
    XSigmaX <- sapply(seq_along(z), function(i) sum(diag(crossprod(z[[i]]) %*% random$Sigma[[i]])))
    wSum <- tapply(w, which, sum)
    sigma2 <- (sum(W*resid(fit)^2) + sum(wSum*XSigmaX))/sum(W)
    wSigma <- add(lapply(seq_along(z), function(i) wSum[i]*random$Sigma[[i]]))
    bb <- add(lapply(seq_along(which), function(i) tcrossprod(random$beta[[i]])*w[i]))
    psi <- (wSigma + bb)/sum(w)
    list(coefficients = if (is(fit, "smooth.spline")) fit$fit else coef(fit),
         sigma2 = list(Random = psi,
           Residual = sigma2),
         df = if (is(fit, "smooth.spline")) fit$df else ncol(x[[1]]))
  }

  object@defineComponent <- expression({
    predict <- function(x, ...) 
      if (is(coef, "smooth.spline.fit")) lapply(x, function(X) stats::predict(coef, X)$y)
      else lapply(x, function(X) X %*% coef)
    
    logLik <- function(x, y, z, which, group, ...) {
      V <- lapply(z, function(Z) tcrossprod(tcrossprod(Z, sigma2$Random), Z) + diag(nrow(Z)) * sigma2$Residual)
      mu <- predict(x, ...)
      llh <- sapply(seq_along(x), function(i) 
                    mvtnorm::dmvnorm(t(y[[i]]), mean = mu[[i]], sigma = V[[which[i]]], log=TRUE)/nrow(V[[which[i]]]))
      as.vector(ungroupPriors(matrix(llh), group, !duplicated(group)))
    }
    new("FLXcomponentlmm",
        parameters = list(coef = coef, sigma2 = sigma2),
        random = list(),
        logLik = logLik, predict = predict,
        df = df)
  })

  determineRandom <- function(mu, y, z, which, sigma2) {
    Sigma <- lapply(z, function(Z)
                    solve(crossprod(Z) / sigma2$Residual + solve(sigma2$Random)))
    Sigma_tilde <-  lapply(seq_along(z), function(i) (tcrossprod(Sigma[[i]], z[[i]])/sigma2$Residual))
    beta <- lapply(seq_along(which), function(i) Sigma_tilde[[which[i]]] %*% (y[[i]] - mu[[i]]))
    list(beta = beta, Sigma = Sigma)
  }
  
  object@fit <- if (any(varFix)) {
    function(x, y, w, z, which, random) {
      fit <- lapply(seq_len(ncol(w)), function(k) lmm.wfit(x, y, w[,k], z, which, random[[k]]))
      if (varFix["Random"]) {
        prior_w <- apply(w, 2, weighted.mean, w = sapply(x, length))
        Random <- add(lapply(seq_along(fit), function(i) fit[[i]]$sigma2$Random * prior_w[i]))
        for (i in seq_along(fit)) fit[[i]]$sigma2$Random <- Random
      }
      if (varFix["Residual"]) {
        prior <- colMeans(w)
        Residual <- sum(sapply(fit, function(x) x$sigma2$Residual) * prior)
        for (i in seq_along(fit)) fit[[i]]$sigma2$Residual <- Residual
      }
      n <- nrow(fit[[1]]$sigma2$Random)
      lapply(fit, function(Z) {
        comp <- with(list(coef = coef(Z),
                          sigma2 =  Z$sigma2,
                          df = Z$df + n*(n+1)/(2*ifelse(varFix["Random"], ncol(w), 1)) + ifelse(varFix["Residual"], 1/ncol(w), 1)),
                     eval(object@defineComponent))
        comp@random <- determineRandom(comp@predict(x), y, z, which, comp@parameters$sigma2)
        comp
      })
    }
  } else {
    function(x, y, w, z, which, random){
      fit <- lmm.wfit(x, y, w, z, which, random)
      n <- nrow(fit$sigma2$Random)
      comp <- with(list(coef = coef(fit),
                        df = fit$df + n*(n+1)/2 + 1,
                        sigma2 =  fit$sigma2),
                   eval(object@defineComponent))
      comp@random <- determineRandom(comp@predict(x), y, z, which, comp@parameters$sigma2)
      comp
    }
  }
  object
}

setMethod("FLXmstep", signature(model = "FLXMRlmm"),
          function(model, weights, components)
{
  weights <- weights[!duplicated(model@group),,drop=FALSE]
  if (!is(components[[1]], "FLXcomponentlmm")) {
    random <- list(beta = lapply(model@which, function(i) rep(0, ncol(model@z[[i]]))),
                   Sigma = lapply(model@z, function(x) diag(ncol(x))))
    return(sapply(seq_len(ncol(weights)),
                  function(k) model@fit(model@x, model@y, weights[,k], model@z, model@which, random)))
 }else {
   return(sapply(seq_len(ncol(weights)),
                 function(k) model@fit(model@x, model@y, weights[,k], model@z, model@which, 
                                       components[[k]]@random)))
 }
})

setMethod("FLXmstep", signature(model = "FLXMRlmmfix"),
          function(model, weights, components)
{
  weights <- weights[!duplicated(model@group),,drop=FALSE]
  if (!is(components[[1]], "FLXcomponentlmm")) {
    random <- rep(list(list(beta = lapply(model@which, function(i) rep(0, ncol(model@z[[i]]))),
                            Sigma = lapply(model@z, function(x) diag(ncol(x))))), ncol(weights))
    return(model@fit(model@x, model@y, weights, model@z, model@which, random))
  }else
   return(model@fit(model@x, model@y, weights, model@z, model@which, lapply(components, function(x) x@random)))
})


setMethod("FLXgetModelmatrix", signature(model="FLXMRlmm"),
          function(model, data, formula, lhs=TRUE, ...)
{
  formula_nogrouping <- RemoveGrouping(formula)
  if (identical(paste(deparse(formula_nogrouping), collapse = ""), paste(deparse(formula), collapse = ""))) stop("please specify a grouping variable")
  model <- callNextMethod(model, data, formula, lhs)
  model@fullformula <- update(model@fullformula,
                              paste(".~. |", .FLXgetGroupingVar(formula)))
  mt1 <- terms(model@random, data=data)
  mf <- model.frame(delete.response(mt1), data=data, na.action = NULL)
  model@z <- model.matrix(attr(mf, "terms"), data)
  model@group <- grouping <- .FLXgetGrouping(formula, data)$group
  model@x <- matrix(lapply(unique(grouping), function(g) model@x[grouping == g, , drop = FALSE]), ncol = 1)
  if (lhs) model@y <- matrix(lapply(unique(grouping), function(g) model@y[grouping == g, , drop = FALSE]), ncol = 1)
  z <- lapply(unique(grouping), function(g) model@z[grouping == g, , drop = FALSE])
  z1 <- unique(z)
  model@which <- sapply(z, function(y) which(sapply(z1, function(x) isTRUE(all.equal(x, y)))))
  model@z <- matrix(z1, ncol = 1)
  model
})

setMethod("FLXgetObs", "FLXMRlmm", function(model) sum(sapply(model@x, nrow)))

setMethod("FLXdeterminePostunscaled", signature(model = "FLXMRlmm"), function(model, components, ...) {
  sapply(components, function(x) x@logLik(model@x, model@y, model@z, model@which, model@group))
})

setMethod("predict", signature(object="FLXMRlmm"), function(object, newdata, components, ...)
{
  object <- FLXgetModelmatrix(object, newdata, formula = object@fullformula, lhs = FALSE)
  lapply(components, function(comp) unlist(comp@predict(object@x, ...)))
})

setMethod("rFLXM", signature(model="FLXMRlmm", components="list"),
          function(model, components, class, group, ...) {
            class <- class[!duplicated(group)]
            y <- NULL
            for (l in seq_along(components)) {
              yl <- as.matrix(rFLXM(model, components[[l]], ...))
              if (is.null(y))  y <- matrix(NA, nrow = length(class), ncol = ncol(yl))
              y[class == l,] <- yl[class==l,,drop=FALSE]
              y <- matrix(y, ncol = ncol(yl))
            }
            y 
          })

setMethod("rFLXM", signature(model = "FLXMRlmm", components = "FLXcomponent"),
          function(model, components, ...) {
            sigma2 <- components@parameters$sigma2
            V <- lapply(model@z, function(Z) tcrossprod(tcrossprod(Z, sigma2$Random), Z) + diag(nrow(Z)) * sigma2$Residual)
            mu <- components@predict(model@x)
            matrix(lapply(seq_along(model@x), function(i) 
                          t(mvtnorm::rmvnorm(1, mean = mu[[i]], sigma = V[[model@which[i]]]))), ncol = 1)
                 })

setMethod("FLXgetNewModelmatrix", "FLXMRlmm", function(object, model, indices, groups) {
  object@y <- model@y[indices,,drop=FALSE]
  object@x <- model@x[indices,,drop=FALSE]
  object@which <- model@which[indices]
  if (length(unique(object@which)) < length(object@z)) {
    object@z <- model@z[sort(unique(object@which)),,drop=FALSE]
    object@which <- match(object@which, sort(unique(object@which)))
  }
  object
})

