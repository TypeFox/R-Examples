setGeneric("boot", function(object, ...) standardGeneric("boot"))
setGeneric("LR_test", function(object, ...) standardGeneric("LR_test"))

setClass("FLXboot",
         representation(call="call",
                        object="flexmix",
                        parameters="list",
                        concomitant="list",
                        priors="list",
                        logLik="matrix",
                        k="matrix",
                        converged="matrix",
                        models="list",
                        weights="list"))

setMethod("show", "FLXboot",
function(object) {
  cat("\nCall:", deparse(object@call,0.75*getOption("width")),
      sep="\n")
})       

generate_weights <- function(object) {
  if(is.null(object@weights) & is(object@model, "FLXMC")) {
    X <- do.call("cbind", lapply(object@model, function(z) z@y))
    x <- apply(X, 1, paste, collapse = "")
    x <- as.integer(factor(x, unique(x)))
    object@weights <- as.vector(table(x))
    indices_unique <- !duplicated(x)
    for (i in seq_along(object@model)) {
      object@model[[i]]@x <- object@model[[i]]@x[indices_unique,,drop=FALSE]
      object@model[[i]]@y <- object@model[[i]]@y[indices_unique,,drop=FALSE]
    }
    object@concomitant@x <- object@concomitant@x[indices_unique,,drop=FALSE]
  }
  object
}

setGeneric("FLXgetNewModelmatrix", function(object, ...) standardGeneric("FLXgetNewModelmatrix"))

setMethod("FLXgetNewModelmatrix", "FLXM", function(object, model, indices, groups) {
  if (length(groups$group) > 0) {
    obs_groups <- lapply(groups$group[groups$groupfirst][indices],
                         function(x) which(x == groups$group))
    indices_grouped <- unlist(obs_groups)
  } else {
    indices_grouped <- indices
  }
  object@y <- model@y[indices_grouped,,drop=FALSE]
  object@x <- model@x[indices_grouped,,drop=FALSE]
  object
})

setMethod("FLXgetNewModelmatrix", "FLXMRglmfix", function(object, model, indices, groups) {
  if (length(groups$group) > 0) {
    obs_groups <- lapply(groups$group[groups$groupfirst][indices],
                         function(x) which(x == groups$group))
    indices_grouped <- unlist(obs_groups)
  } else {
    indices_grouped <- indices
  }
  object@y <- do.call("rbind", rep(list(model@y[indices_grouped,,drop=FALSE]), sum(model@nestedformula@k)))
  object@x <- do.call("rbind", lapply(seq_len(sum(model@nestedformula@k)), function(K)
                                      model@x[model@segment[,K],,drop=FALSE][indices_grouped,,drop=FALSE]))
  N <- nrow(object@x)/sum(model@nestedformula@k)
  object@segment <- matrix(FALSE, ncol = sum(model@nestedformula@k), nrow = nrow(object@x))
  for (m in seq_len(sum(model@nestedformula@k))) object@segment[(m - 1) * N + seq_len(N), m] <- TRUE
  object
})

boot_flexmix <- function(object, R, sim = c("ordinary", "empirical", "parametric"), initialize_solution = FALSE,
         keep_weights = FALSE, keep_groups = TRUE, verbose = 0, control, k, model = FALSE, ...) {
  sim <- match.arg(sim)
  if (missing(R)) stop("R needs to be specified")
  if (!missing(control)) object@control <- do.call("new", c(list(Class = "FLXcontrol", object@control), control))
  if (missing(k)) k <- object@k
  m <- length(object@model)
  has_weights <- !keep_weights & !is.null(object@weights)
  if (has_weights) object <- undo_weights(object)
  if (!keep_groups & length(object@group)) {
    object@concomitant@x <- object@concomitant@x[as.integer(object@group),,drop = FALSE]
    object@group <- factor()
  }
  groups <- list()
  groups$group <- object@group
  groups$groupfirst <- if (length(groups$group) > 0) groupFirst(groups$group)
                       else rep(TRUE, FLXgetObs(object@model[[1]]))

  concomitant <- parameters <- priors <- models <- weights <- vector("list", R)
  logLik <- ks <- converged <- matrix(nrow=R, ncol = length(k), dimnames = list(BS = seq_len(R), k = k))
  for (iter in seq_len(R)) {
    new <- object
    newgroups <- groups
    if(verbose && !(iter%%verbose)) cat("* ")
    if (iter > 1) {
      if (sim == "parametric") {
        y <- rflexmix(object, ...)$y
        for (i in seq_len(m))
          new@model[[i]]@y <- matrix(as.vector(t(y[[i]])),
                                     nrow = nrow(new@model[[i]]@x),
                                     ncol = ncol(y[[i]]), byrow = TRUE)
      } else {
        n <- sum(groups$groupfirst)
        indices <- sample(seq_len(n), n, replace = TRUE)
        if (length(groups$group) > 0) {
            obs_groups <- lapply(groups$group[groups$groupfirst][indices],
                                 function(x) which(x == groups$group))
            newgroups$group <- factor(rep(seq_along(obs_groups), sapply(obs_groups, length)))
            newgroups$groupfirst <- !duplicated(newgroups$group)
        }
        for (i in seq_len(m)) {
          new@model[[i]] <- FLXgetNewModelmatrix(new@model[[i]], object@model[[i]],
                                                 indices, groups)
        }
        new@concomitant@x <- new@concomitant@x[indices,,drop=FALSE]
      }
    }
    if (has_weights & !length(groups$group) > 0) {
      new <- generate_weights(new)
      newgroups$groupfirst <- rep(TRUE, FLXgetObs(new@model[[1]]))
    }
    parameters[[iter]] <- concomitant[[iter]] <- priors[[iter]] <- list()
    NREP <- rep(object@control@nrep, length(k))
    if (initialize_solution & object@k %in% k) NREP[k == object@k] <- 1L
    for (K in seq_along(k)) {
      fit <- new("flexmix", logLik = -Inf)
      for (nrep in seq_len(NREP[K])) {
        if (k[K] != object@k | !initialize_solution)  {
          postunscaled <- initPosteriors(k[K], NULL, FLXgetObs(new@model[[1]]), newgroups)
        } else {
          postunscaled <- matrix(0, nrow = FLXgetObs(new@model[[1]]), ncol = k[K])
          for (i in seq_len(m)) 
            postunscaled <- postunscaled + FLXdeterminePostunscaled(new@model[[i]], lapply(new@components, function(x) x[[i]]))
          if(length(newgroups$group)>0)
            postunscaled <- groupPosteriors(postunscaled, newgroups$group)
          prior <- evalPrior(new@prior, new@concomitant)
          postunscaled <- if (is(prior, "matrix")) postunscaled + log(prior)
                          else sweep(postunscaled, 2, log(prior), "+")
          postunscaled <- exp(postunscaled - log_row_sums(postunscaled))
        }
        x <- try(FLXfit(new@model, new@concomitant, new@control, postunscaled, newgroups, weights = new@weights))
        if (!is(x, "try-error")) {
          if(logLik(x) > logLik(fit))
            fit <- x
        }
      }
      if (is.finite(logLik(fit))) {
        parameters[[iter]][paste(k[K])] <- list(parameters(fit, simplify = FALSE, drop = FALSE))
        concomitant[[iter]][paste(k[K])] <- list(parameters(fit, which = "concomitant"))
        priors[[iter]][[paste(k[K])]] <- prior(fit)
        logLik[iter, paste(k[K])] <- logLik(fit)
        ks[iter, paste(k[K])] <- fit@k
        converged[iter, paste(k[K])] <- fit@converged
        if (model) {
          models[[iter]] <- fit@model
          weights[[iter]] <- fit@weights
        }
      } else {
        parameters[[iter]][[paste(k[K])]] <- concomitant[[iter]][[paste(k[K])]] <- priors[[iter]][[paste(k[K])]] <- NULL
      }
    }
  }
  if(verbose) cat("\n")
  new("FLXboot", call = sys.call(-1), object = object, parameters = parameters,
      concomitant = concomitant, priors = priors, logLik = logLik, k = ks,
      converged = converged, models = models, weights = weights)
}

setMethod("boot", signature(object="flexmix"), boot_flexmix)

setMethod("LR_test",
          signature(object="flexmix"),
function(object, R, alternative = c("greater", "less"), control, ...) {
  alternative <- match.arg(alternative)
  if (missing(control)) control <- object@control
  if (object@k == 1 & alternative == "less") stop(paste("alternative", alternative, "only possible for a mixture\n",
                        "with at least two components"))
  k <- object@k + switch(alternative, greater = 0:1, less = 0:-1)
  names(k) <- k
  boot <- boot(object, R, sim = "parametric", k = k, 
                    initialize_solution = TRUE, control = control, ...)
  ok <- apply(boot@k, 1, identical, k)
  lrts <- 2*apply(boot@logLik[ok,order(k)], 1, diff)
  STATISTIC <- lrts[1]
  names(STATISTIC) <- "LRTS"
  PARAMETER <- length(lrts)
  names(PARAMETER) <- "BS"
  RETURN <- list(parameter = PARAMETER,
                 p.value = sum(lrts[1] <= lrts)/length(lrts),
                 alternative = alternative,
                 null.value = object@k,
                 method = "Bootstrap likelihood ratio test",
                 data.name = deparse(substitute(object)),
                 bootstrap.results = boot)
  class(RETURN) <- "htest"
  RETURN
})

setMethod("parameters", "FLXboot", function(object, k, ...) {
  if (missing(k)) k <- object@object@k
  Coefs <- lapply(seq_along(object@parameters), function(i) 
                  if (is.na(object@k[i])) NULL
                  else do.call("cbind", c(lapply(seq_len(object@k[i]), function(j) 
                                                 unlist(sapply(seq_along(object@object@model), function(m) 
                                                               FLXgetParameters(as(object@object@model[[m]], "FLXMR"),
                                                                                list(with(c(object@parameters[[i]][[paste(k)]][[m]][[j]],
                                                                                            list(df = object@object@components[[j]][[m]]@df)),
                                                                                          eval(object@object@model[[m]]@defineComponent))))))),
                                          as.list(rep(NA, k - object@k[i])))))
  Coefs <- t(do.call("cbind", Coefs))
  colnames(Coefs) <- gsub("Comp.1_", "", colnames(Coefs))
  Prior <- t(do.call("cbind", lapply(object@concomitant,
                                     function(x) do.call("cbind", c(list(x[[paste(k)]]),
                                                                    as.list(rep(NA, k - ifelse(length(x), ncol(x[[paste(k)]]), k))))))))
  cbind(Coefs, Prior)
})

setMethod("clusters", signature(object = "FLXboot", newdata = "listOrdata.frame"), function(object, newdata, k, ...) {
  if (missing(k)) k <- object@object@k
  lapply(seq_len(length(object@priors)), function(i) {
    new <- object@object
    new@prior <- object@priors[[i]][[paste(k)]]
    new@k <- length(new@prior)
    new@components <- rep(list(vector("list", length(object@object@model))), length(new@prior))
    for (m in seq_along(new@model)) {
      variables <- c("x", "y", "offset", "family")
      variables <- variables[variables %in% slotNames(new@model[[m]])]
      for (var in variables) assign(var, slot(new@model[[m]], var))
      for (K in seq_len(object@k[i])) {
        new@components[[K]][[m]] <- with(c(object@parameters[[i]][[paste(k)]][[m]][[K]],
                                           list(df = object@object@components[[K]][[m]]@df)),
                                         eval(object@object@model[[m]]@defineComponent))
      }
    }
    clusters(new, newdata = newdata)})
})

setMethod("posterior", signature(object = "FLXboot", newdata = "listOrdata.frame"), function(object, newdata, k, ...) {
  if (missing(k)) k <- object@object@k
  lapply(seq_len(length(object@priors)), function(i) {
    new <- object@object
    new@prior <- object@priors[[i]][[paste(k)]]
    new@k <- length(new@prior)
    new@components <- rep(list(vector("list", length(object@object@model))), length(new@prior))
    for (m in seq_along(new@model)) {
      variables <- c("x", "y", "offset", "family")
      variables <- variables[variables %in% slotNames(new@model[[m]])]
      for (var in variables) assign(var, slot(new@model[[m]], var))
      for (K in seq_len(object@k[i])) {
        new@components[[K]][[m]] <- with(c(object@parameters[[i]][[paste(k)]][[m]][[K]],
                                           list(df = object@object@components[[K]][[m]]@df)),
                                         eval(object@object@model[[m]]@defineComponent))
      }
    }
    posterior(new, newdata = newdata)})
})

setMethod("predict", signature(object = "FLXboot"), function(object, newdata, k, ...) {
  if (missing(k)) k <- object@object@k
  lapply(seq_len(length(object@priors)), function(i) {
    new <- object@object
    new@components <- vector("list", object@k[i, paste(k)])
    new@components <- lapply(new@components, function(x) vector("list", length(new@model)))
    for (m in seq_along(new@model)) {
      variables <- c("x", "y", "offset", "family")
      variables <- variables[variables %in% slotNames(new@model[[m]])]
      for (var in variables) assign(var, slot(new@model[[m]], var))
      for (K in seq_len(object@k[i, paste(k)])) {
        new@components[[K]][[m]] <- with(c(object@parameters[[i]][[paste(k)]][[m]][[K]],
                                           list(df = object@object@components[[1]][[m]]@df)),
                                         eval(object@object@model[[m]]@defineComponent))
      }
    }
    predict(new, newdata = newdata, ...)})
})

