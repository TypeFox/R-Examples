setClass("FLXMRlmer",
         representation(random = "formula", 
                        lmod = "list",
                        control = "ANY",
                        preproc.z = "function"),
         prototype(preproc.z = function(x, ...) x),
         contains = "FLXMRglm")

defineComponent_lmer <- expression({
    predict <- function(x, ...) x%*%coef
    logLik <- function(x, y, lmod, ...) {
      z <- as.matrix(lmod$reTrms$Zt)
      grouping <- lmod$reTrms$flist[[1]]
      llh <- vector(length=nrow(x))
      for (i in seq_len(nlevels(grouping))) {
        index1 <- which(grouping == levels(grouping)[i])
        index2 <- rownames(z) %in% levels(grouping)[i]
        V <- crossprod(z[index2,index1,drop=FALSE], sigma2$Random) %*% z[index2, index1, drop=FALSE] + diag(length(index1)) * sigma2$Residual
        llh[index1] <- mvtnorm::dmvnorm(y[index1,], mean=predict(x[index1,,drop=FALSE], ...), sigma = V, log=TRUE)/length(index1)
      }
      llh
    }
      
    new("FLXcomponent",
        parameters=list(coef=coef, sigma2=sigma2),
        logLik=logLik, predict=predict,
        df=df)
  })

FLXMRlmer <- function(formula = . ~ ., random, weighted = TRUE, 
                      control = list(), eps = .Machine$double.eps)
{
  random <- if (length(random) == 3) random else formula(paste(".", paste(deparse(random), collapse = "")))
  missCtrl <- missing(control)
  if (missCtrl || !inherits(control, "lmerControl")) {
      if (!is.list(control)) 
          stop("'control' is not a list; use lmerControl()")
      control <- do.call(lme4::lmerControl, control)
  }
  
  object <- new("FLXMRlmer", formula = formula, random = random, control = control,
                family = "gaussian", weighted = weighted, name = "FLXMRlmer:gaussian")
  if (weighted) object@preproc.z <- function(lmod) { 
    if (length(unique(names(lmod[["reTrms"]][["flist"]]))) != 1) stop("only a single variable for random effects is allowed")
    for (i in seq_along(lmod[["reTrms"]][["flist"]])) {
      DIFF <- t(sapply(levels(lmod[["reTrms"]]$flist[[i]]), function(id) {
        index1 <- which(lmod[["reTrms"]]$flist[[i]] == id)
        index2 <- rownames(lmod[["reTrms"]]$Zt) == id
        sort(apply(lmod[["reTrms"]]$Zt[index2, index1, drop=FALSE], 1, paste, collapse = ""))
      }))
      if (length(unique(table(lmod[["reTrms"]][["flist"]][[i]]))) != 1 || nrow(unique(DIFF)) != 1)
        stop("FLXMRlmer does only work correctly if the covariates of the random effects are the same for all observations")
    }
    lmod
  }

  lmer.wfit <- function(x, y, w, lmod) {
    zero.weights <- any(w < eps)
    if (zero.weights) {
      ok <- w >= eps
      w <- w[ok]
      lmod[["fr"]] <- lmod[["fr"]][ok,,drop=FALSE]
      lmod[["X"]] <- lmod[["X"]][ok,,drop=FALSE]
      lmod[["reTrms"]][["Zt"]] <- lmod[["reTrms"]][["Zt"]][,ok]
      lmod[["reTrms"]][["flist"]] <- lmod[["reTrms"]][["flist"]][ok,,drop=FALSE]
    }
    wts <- sqrt(w)
    lmod$X <- lmod$X * wts
    lmod$fr[[1]] <- lmod$fr[[1]] * wts
    devfun <- do.call(lme4::mkLmerDevfun, c(lmod, list(start = NULL, verbose = FALSE, 
                                                 control = control)))
    opt <- lme4::optimizeLmer(devfun, optimizer = control$optimizer, 
                              restart_edge = control$restart_edge, control = control$optCtrl, 
                              verbose = FALSE, start = NULL)
    mer <- lme4::mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
    sigma_res <- lme4::sigma(mer) / sqrt(mean(w))
    vc <- lme4::VarCorr(mer)
    n <- c(0, cumsum(sapply(vc, ncol)))
    Random <- matrix(0, max(n), max(n))
    for (i in seq_along(vc)) {
        index <- (n[i]+1):n[i+1]
        Random[index, index] <- vc[[i]]
    }
    Random <- Random / mean(w)
    list(coefficients = lme4::fixef(mer),
         sigma2 = list(Random  = Random,
             Residual = sigma_res^2),
         df = length(lme4::fixef(mer)) + 1 + length(mer@theta))
  }
  
  object@defineComponent <- defineComponent_lmer
  
  object@fit <- function(x, y, w, lmod){
    fit <- lmer.wfit(x, y, w, lmod)
    with(list(coef = coef(fit),
              df = fit$df,
              sigma2 =  fit$sigma2),
         eval(object@defineComponent))
  }
  object
}

setMethod("FLXgetModelmatrix", signature(model="FLXMRlmer"),
          function(model, data, formula, lhs=TRUE, contrasts = NULL, ...)
{
  formula_nogrouping <- RemoveGrouping(formula)
  if (identical(paste(deparse(formula_nogrouping), collapse = ""), paste(deparse(formula), collapse = ""))) stop("please specify a grouping variable")
  model <- callNextMethod(model, data, formula, lhs)
  random_formula <- update(model@random,
                           paste(".~. |", .FLXgetGroupingVar(formula)))
  fullformula <- model@fullformula
  if (!lhs) fullformula <- fullformula[c(1,3)]
  fullformula <- update(fullformula,
                        paste(ifelse(lhs, ".", ""), "~. + ", paste(deparse(random_formula[[3]]), collapse = "")))
  model@fullformula <- update(model@fullformula,
                              paste(ifelse(lhs, ".", ""), "~. |", .FLXgetGroupingVar(formula)))
  model@lmod <- lme4::lFormula(fullformula, data, REML = FALSE, control = model@control)
  model@lmod <- model@preproc.z(model@lmod)
  model
})

setMethod("FLXmstep", signature(model = "FLXMRlmer"),
          function(model, weights, ...)
{
  apply(weights, 2, function(w) model@fit(model@x, model@y, w, model@lmod))
})

setMethod("FLXdeterminePostunscaled", signature(model = "FLXMRlmer"), function(model, components, ...) {
  sapply(components, function(x) x@logLik(model@x, model@y, model@lmod))
})

