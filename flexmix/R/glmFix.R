#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: glmFix.R 4978 2014-02-13 15:45:15Z gruen $
#

FLXMRglmfix <- function(formula=.~., fixed=~0, varFix = FALSE, nested = NULL,
                        family=c("gaussian", "binomial", "poisson", "Gamma"),
                        offset=NULL)
{
      family <- match.arg(family)
      nested <- as(nested, "FLXnested")
      
      if (length(fixed) == 3) stop("no left hand side allowed for fixed")       
      z <- new("FLXMRglmfix", FLXMRglm(formula, family, offset),
               fixed=fixed, name=paste("FLXMRglmfix", family, sep=":"), nestedformula=nested,
               variance = varFix)

      if(family=="gaussian"){
        z@fit <- function(x, y, w, incidence, variance, ...){
            fit <- lm.wfit(x, y, w=w, offset=offset)
            k <- nrow(incidence)
            n <- nrow(x)/k
            sigma <- vector(length=k)
            cumVar <- cumsum(c(0, variance))
            for (i in seq_along(variance)) {
              ind <- cumVar[i]*n + seq_len(n*variance[i])
              sigma[cumVar[i] + seq_len(variance[i])] <- sqrt(sum(fit$weights[ind] * fit$residuals[ind]^2 /
                                                           mean(fit$weights[ind]))/ (length(ind) - sum(incidence[i,])))
            }
            fit <- fit[c("coefficients")]
            coefs <- coef(fit)
            names(coefs) <- colnames(incidence)
            df <- rowSums(incidence/rep(colSums(incidence), each = nrow(incidence))) + rep(1/variance, variance)
            lapply(seq_len(k),
                   function(K) with(list(coef=coefs[as.logical(incidence[K,])],
                                         sigma=sigma[K],
                                         df= df[K]),
                                    eval(z@defineComponent)))
          }
    }
    else if(family=="binomial"){
        z@fit <- function(x, y, w, incidence, ...){
            fit <- glm.fit(x, y, weights=w, family=binomial(), offset=offset)
            fit <- fit[c("coefficients","family")]
            k <- nrow(incidence)
            coefs <- coef(fit)
            names(coefs) <- colnames(incidence)
            df <- rowSums(incidence/rep(colSums(incidence), each = nrow(incidence)))
            lapply(seq_len(k),
                   function(K) with(list(coef=coefs[as.logical(incidence[K,])],
                                         df = df[K]),
                                    eval(z@defineComponent)))
          }
    }
    else if(family=="poisson"){
        z@fit <- function(x, y, w, incidence, ...){
            fit <- glm.fit(x, y, weights=w, family=poisson(), offset=offset)
            fit <- fit[c("coefficients","family")]
            k <- nrow(incidence)
            coefs <- coef(fit)
            names(coefs) <- colnames(incidence)
            df <- rowSums(incidence/rep(colSums(incidence), each = nrow(incidence)))
            lapply(seq_len(k),
                   function(K) with(list(coef=coefs[as.logical(incidence[K,])],
                                         df = df[K]),
                                    eval(z@defineComponent)))
          }
    }
    else if(family=="Gamma"){
        z@fit <- function(x, y, w, incidence, ...){
            fit <- glm.fit(x, y, weights=w, family=Gamma(), offset=offset)
            shape <- sum(fit$prior.weights)/fit$deviance
            fit <- fit[c("coefficients","family")]
            k <- nrow(incidence)
            coefs <- coef(fit)
            names(coefs) <- colnames(incidence)
            df <- rowSums(incidence/rep(colSums(incidence), each = nrow(incidence)))
            lapply(seq_len(k),
                   function(K) with(list(coef=coefs[as.logical(incidence[K,])],
                                         df = df[K],
                                         shape = shape),
                                    eval(z@defineComponent)))
          }
    }
    else stop(paste("Unknown family", family))
    z
}

###**********************************************************

setMethod("refit_mstep", signature(object="FLXMRglmfix", newdata="missing"),
function(object, newdata, weights, ...)
{
    z <- new("FLXRMRglmfix", design=object@design)
    z@fitted <- object@refit(object@x,
                             object@y,
                             as.vector(weights))
    z <- rep(list(z), nrow(object@design))
    for (k in seq_len(nrow(object@design))) z[[k]]@design <- object@design[k,, drop=FALSE]
    z
})

###**********************************************************

setMethod("fitted", signature(object="FLXMRglmfix"),
function(object, components, ...)
{
    N <- nrow(object@x)/length(components)
    z <- list()
    for(n in seq_along(components)){
      x <- object@x[(n-1)*N + seq_len(N), as.logical(object@design[n,]), drop=FALSE]
      z[[n]] <- list(components[[n]]@predict(x))
    }
    z
})
               
###**********************************************************

setMethod("predict", signature(object="FLXMRglmfix"),
function(object, newdata, components, ...)
{
  model <- FLXgetModelmatrix(object, newdata, object@fullformula, lhs=FALSE)
  k <- sum(object@nestedformula@k)
  N <- nrow(model@x)/k
  z <- list()
  for (m in seq_len(k)) {
    z[[m]] <- components[[m]]@predict(model@x[model@segment[,m], as.logical(model@design[m,]), drop=FALSE], ...)
  }
  z
})

###**********************************************************

setMethod("FLXgetModelmatrix", signature(model="FLXMRfix"), function(model, data, formula, lhs=TRUE, ...)
{
  formula <- RemoveGrouping(formula)
  if (length(grep("\\|", deparse(model@formula)))) stop("no grouping variable allowed in the model")
    if(is.null(model@formula))
      model@formula <- formula
    model@fullformula <- update.formula(formula, model@formula)
    k <- model@nestedformula
    mm.all <- modelMatrix(model@fullformula, model@fixed, k@formula, data, lhs, model@xlevels)
    model@design <- modelDesign(mm.all, k)
    desNested <- if (sum(sapply(mm.all$nested, ncol))) {
      rbind(ncol(mm.all$fixed) + seq_len(sum(sapply(mm.all$nested, ncol))),
            unlist(lapply(seq_along(mm.all$nested), function(i) rep(i, ncol(mm.all$nested[[i]])))))
    }else  matrix(ncol=0, nrow=2)
    model@x <- cbind(kronecker(rep(1, sum(k@k)), mm.all$fixed),
                     do.call("cbind", lapply(unique(desNested[2,]), function(i) {
                       kronecker(model@design[,desNested[1, desNested[2, ] == i][1]],
                                 mm.all$nested[[i]])})),
                     kronecker(diag(sum(k@k)), mm.all$random))
    N <- nrow(model@x)/sum(k@k)
    model@segment <- matrix(FALSE, ncol = sum(k@k), nrow = nrow(model@x))
    for (m in seq_len(sum(k@k))) model@segment[(m - 1) * N + seq_len(N), m] <- TRUE
    if (lhs) {
      y <- mm.all$response
      rownames(y) <- NULL
      response <- as.matrix(apply(y, 2, rep, sum(k@k)))
      model@y <- model@preproc.y(response)
    }
    model@x <- model@preproc.x(model@x)
    model@xlevels <- mm.all$xlevels
    model
})

