## interface 
btmix <- function(formula, data, k, subset, weights,
                     nrep = 3, cluster = NULL, control = NULL,
                     verbose = TRUE, drop = TRUE, unique = FALSE, which = NULL,
		     type = c("loglin", "logit"), ref = NULL, undecided = NULL, position = NULL, ...){  
  ## process call
  cl <- match.call()
  has_subset <- !missing(subset)
  has_weights <- !missing(weights)

  ## arrange formula and data arguments correctly
  if(missing(formula)) {
    d <- as.data.frame(matrix(0, nrow = nrow(data), ncol = 0L))
    d$.response <- as.matrix(data)
    conc <- NULL 
    pc <- data
  } else if(!inherits(formula, "formula") & missing(data)) {
    d <- as.data.frame(matrix(0, nrow = nrow(formula), ncol = 0L))
    d$.response <- as.matrix(formula)
    conc <- NULL
    pc <- formula
  } else {
  
    ## process call
    if(missing(data)) data <- environment(formula)
    aux <- match.call(expand.dots = FALSE)
    aux <- aux[c(1L, match(c("formula", "data", "subset", "weights"), names(aux), 0L))]
    
    ## process formula via Formula
    ff <- Formula(formula)

    ## handle conomitants first
    auxc <- aux
    auxc[[1]] <- as.name("get_all_vars")
    auxc$formula <- formula(ff, lhs = 0, rhs = NULL)
    if(has_subset) names(auxc)[names(auxc) == "subset"] <- ".subset"
    if(has_weights) names(auxc)[names(auxc) == "weights"] <- ".weights"
    d <- eval(auxc, parent.frame())
    if(has_subset) d <- d[d$.subset, , drop = FALSE]
    d$.subset <- NULL
    conc <- ncol(d) > has_weights
    conc <- if(conc) {
      f <- formula(ff, lhs = 0, rhs = NULL, collapse = TRUE)
      FLXPmultinom(f)
    } else {
      NULL
    }
    
    ## handle response
    auxr <- aux
    auxr[[1]] <- as.name("model.frame")
    auxr$formula <- ff
    auxr$lhs <- NULL
    auxr$rhs <- 0
    auxr$na.action <- na.pass
    auxr$weights <- NULL
    mf <- eval(auxr, parent.frame())
    pc <- model.part(ff, lhs = NULL, rhs = 0, data = mf)[[1L]]
    d$.response <- as.matrix(pc)
  }
  
  ## FIXME: look at NA processing
  ## FIXME: store paircomp attributes in btmix object
  
  ## data processing: remove observations without any paired comparisons
  ## missing.obs <- is.na(d$.response)
  ## d <- d[!missing.obs, , drop = FALSE]

  ## data processing: remove observations with any NA (required for flexmix)
  # (removing only observations with missings in covariates is not enough)
##  cc <- complete.cases(d)
##  d <- d[cc, , drop = FALSE]
  
  ## driver function
  model <- FLXMCbtreg(type = type, ref = ref, undecided = undecided, position = position)

  ## control parameters
  ctrl <- as(control, "FLXcontrol")  
    
  ## fitting
  z <- if(has_weights) {
    stepFlexmix(.response ~ 1, data = d, k = k, nrep = nrep, cluster = cluster, model = model,
                concomitant = conc, control = ctrl, verbose = verbose,
                drop = drop, unique = unique, weights = d$.weights, ...)
  } else {
    stepFlexmix(.response ~ 1, data = d, k = k, nrep = nrep, cluster = cluster, model = model,
                concomitant = conc, control = ctrl, verbose = verbose,
                drop = drop, unique = unique, ...)
  }

  ## select model (if desired)
  if (!is.null(which) & class(z) == "stepFlexmix"){
    z <- getModel(z, which = which)
  }
  
  ## classify
  if (class(z) == "flexmix") {
    z <- as(z, "btmix")
    z@flx.call <- z@call
    z@call <- cl
    z@labels <- labels(pc)
    z@mscale <- mscale(pc)
    z@undecided <- if(is.null(undecided)) z@mscale[2] == 0 else undecided
    z@ref <- ref
    z@type <- match.arg(type, c("loglin", "logit"))
  }
  else {
    z@models <- lapply(z@models, function(model){
      model <- as(model, "btmix")
      model@flx.call <- model@call
      model@call <- cl
      model@call$k <- model@flx.call$k
      model@labels <- labels(pc)
      model@mscale <- mscale(pc)
      model@undecided <- if(is.null(undecided)) model@mscale[2] == 0 else undecided
      model@ref <- ref
      model@type <- match.arg(type, c("loglin", "logit"))
      return(model)
    })
    z <- as(z, "stepBTmix")
    z@flx.call <- z@call
    z@call <- cl
    z@labels <- labels(pc)
    z@mscale <- mscale(pc)
  }

  return(z)
}



## Flexmix driver for Rasch mixture model
FLXMCbtreg <- function(formula = . ~ ., type = c("loglin", "logit"), ref = NULL,
  undecided = NULL, position = NULL, ...)
{
  retval <- new("FLXMCbt", weighted = TRUE, formula = formula,
                name = "Bradley-Terry mixture model")

  retval@defineComponent <- expression({

    logLik <- function(x,y,...) {
      loglikfun_btreg(coef, y, labels, undecided, ref, type,...)
    }

    new("FLXcomponent", df = df , logLik = logLik,
        parameters = list(coef = coef))
  })

  retval@fit <- function(x,y,w, ...){
    if(!inherits(y, "paircomp")) y <- paircomp(y,
      labels = unique(unlist(strsplit(colnames(y), ":", fixed = TRUE))))

    btreg <- btmodel(y, weights = w, type = type, ref = ref,
      undecided = undecided, position = position, ...)

    para <- list(coef = btreg$coefficients, df = btreg$df,
                 labels = btreg$labels, undecided = btreg$undecided,
                 ref = btreg$ref, type = btreg$type)

    with(para, eval(retval@defineComponent))
  }
  retval@dist <- "Bradley-Terry"
  retval
}

## compute individual contributions to log-likelihood function
loglikfun_btreg <- function(par, data, labels, undecided, ref, type, ...)
{
  nobj <- length(labels)
  npar <- nobj - !undecided
  npc <- nobj * (nobj - 1)/2
  ix <- which(upper.tri(diag(nobj)), arr.ind = TRUE)
  if(is.null(ref)) ref <- nobj
  if(is.character(ref)) ref <- match(ref, labels)
  
  par2logprob <- switch(type,
    "loglin" = function(i) {
      p <- rep(0, nobj)
      p[-ref] <- if(undecided) par[-npar] else par
      p <- p[ix[i,]]
      if(undecided) p <- c(p, par[npar] + mean(p))
      p - log(sum(exp(p)))
    },
    "logit" = function(i) {
      p <- rep(0, nobj)
      p[-ref] <- par
      plogis(c(-1, 1) * diff(p[ix[i,]]), log.p = TRUE)
    }
  )
  logp <- t(sapply(1:npc, par2logprob))

  data01 <- as.matrix(data)
  data01 <- matrix(as.numeric(cbind(
    data01 == 1, data01 == -1, if(undecided) data01 == 0 else NULL)),
    nrow = nrow(data01))
  rval <- drop(data01 %*% as.vector(logp))

  ## return result
  return(rval)
}
