## interface 
raschmix <- function(formula, data, k, subset, weights,
                     scores = c("saturated", "meanvar"), restricted = FALSE,
                     nrep = 3, cluster = NULL, control = NULL,
                     verbose = TRUE, drop = TRUE, unique = FALSE, which = NULL,
		     reltol = 1e-10, deriv = "sum", hessian = FALSE, restart = TRUE,
		     model = NULL, gradtol = reltol, ...){  
  ## process call
  cl <- match.call()
  has_subset <- !missing(subset)
  has_weights <- !missing(weights)
  scores <- match.arg(head(tolower(scores), 1L), c("saturated", "meanvar", "constant"))
  misRes <- missing(restricted)
  
  ## gradtol renamed reltol in RaschModel.fit()/raschmodel()
  if(missing(reltol) && !missing(gradtol) && !is.null(gradtol)) reltol <- gradtol

  ## arrange formula and data arguments correctly
  if(missing(formula)) {
    d <- as.data.frame(matrix(0, nrow = nrow(data), ncol = 0L))
    d$.response <- data
    conc <- NULL 
  } else if(!inherits(formula, "formula") & missing(data)) {
    d <- as.data.frame(matrix(0, nrow = nrow(formula), ncol = 0L))
    d$.response <- formula
    conc <- NULL
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
    #d$.response <- as.matrix(model.part(ff, lhs = NULL, rhs = 0, data = mf))
    mp <- model.part(ff, lhs = NULL, rhs = 0, data = mf)
    if(all(sapply(mp, is.itemresp))) {
      d$.response <- mp[[1]]
      if(length(mp) > 1) for (i in 2:length(mp)) d$.response <- merge(d$.response, mp[[i]])
    } else {
      d$.response <- as.matrix(mp)
    }
  }
  
  ## data processing: remove observations without any item responses
  if(class(d$.response) == "itemresp") {
      missing.obs <- is.na(d$.response)
  } else {
    missing.obs <- apply(is.na(d$.response), 1, all)
  }
  d <- d[!missing.obs, , drop = FALSE]

  ## data processing: remove observations with any NA in concomitant variables
  #conc.na <- apply(is.na(d[-which(".response" == names(d))]), 2, any)
  #d <- d[-conc.na, , drop = FALSE]
  conc.cc <- complete.cases(d[-which(".response" == names(d))])
  d <- d[conc.cc, , drop = FALSE]

  ## data processing: remove observations with any NA (required for flexmix)
  # (removing only observations with missings in covariates is not enough)
  #cc <- complete.cases(d)
  #d <- d[cc, , drop = FALSE]

  ## data processing: non-identified items
  n.total <- nrow(d$.response)
  cm <- colMeans(d$.response, na.rm = TRUE)
  status <- as.character(cut(cm, c(-Inf, 1/(2 * n.total), 1 - 1/(2 * n.total),
                                   Inf), labels = c("0", "0/1", "1")))
  status[is.na(status)] <- "NA"
  status <- factor(status, levels = c("0/1", "0", "1", "NA"))
  if(is.null(colnames(d$.response)))
    colnames(d$.response) <- paste("Item", gsub(" ", "0",
                                  format(1:ncol(d$.response))), sep = "")
  names(status) <- colnames(d$.response)
  ident <- status == "0/1"
  d$.response <- d$.response[,ident, drop = FALSE]
 
  
  ## data processing: extreme scorers
  ## extreme scorer = subjects who score all non-missing items with either 0 or 1
  score.0 <- rowSums(d$.response, na.rm = TRUE) == 0
  n.0 <- sum(score.0)
  d <- d[!score.0, , drop = FALSE]
  if(class(d$.response) == "itemresp"){
    score.m <- (rowSums(d$.response, na.rm = TRUE) + rowSums(is.na(as.matrix(d$.response)))) == ncol(d$.response)
  } else {
    score.m <- (rowSums(d$.response, na.rm = TRUE) + rowSums(is.na(d$.response))) == ncol(d$.response)
  }
  n.m <- sum(score.m)
  d <- d[!score.m, , drop = FALSE]

  pi.0 <- n.0 / n.total
  pi.m <- n.m / n.total
  pi.nonex <- nrow(d$.response) / n.total

  ## reference category for saturated score model
  rs <- factor(rowSums(d$.response, na.rm = TRUE), levels = 1:(ncol(d$.response)-1L))
  rs <- table(rs)
  ref <- min(which(rs > 0))

  ## score model
  if(scores == "constant" & !restricted){
    if(!misRes) warning("Constant score model only available with restricted components.")
    restricted <- TRUE
  }

  ## delta are the parameters of the score model
  if (restricted == TRUE){
    ## calculate score model and df
    scoreMod <- scoreModel(y = d$.response, w = rep(1, nrow(d$.response)), ref = ref, scores = scores)
    delta <- scoreMod$score
    dfScore <- scoreMod$df
  } else {
    delta <- NULL
    dfScore <- 0
  }

  ## Rasch driver including score model
  if(is.null(model)) model <- FLXMCrasch(scores = scores, delta = delta, ref = ref,
     nonExtremeProb = pi.nonex, reltol = reltol, deriv = deriv, hessian = hessian,
     restart = restart)

  ## control parameters
  ctrl <- as(control, "FLXcontrol")  
    
  ## starting values via mrm() from "mRm"
  if(identical(cluster, "mrm")) {
    if(length(k) > 1L) warning("starting values from mrm() can only be used with a single 'k', first used")
    k <- k[1L]
    cluster <- if(nrep > 1L) {
      fits <- lapply(1:nrep, function(i) mRm::mrm(d$.response, k))
      ix <- which.max(sapply(fits, "[[", "logLik"))
      mrm_clusters(d$.response, k, fits[[ix]])
    } else {
      mrm_clusters(d$.response, k)
    }
    nrep <- 1L
  }

  ## fitting
  ## FIXME: use initFlexmix instead of stepFlexmix?
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
  if (class(z) == "flexmix") { ## is(z, "flexmix")?
    ## df: include df for extreme score probs if estimated
    z@df <- z@df + (pi.0 != 0) + (pi.m != 0)
    ## df: include df for score model if restricted over components
    if(restricted) z@df <- z@df + dfScore
    ## logLik: include term for extreme scorers if estimated
    if (pi.0 != 0) z@logLik <- z@logLik + n.0*log(pi.0)
    if (pi.m != 0) z@logLik <- z@logLik + n.m*log(pi.m)
    ## classify as raschmix and fill raschmix-specific slots
    z <- as(z, "raschmix")
    z@scores <- scores
    z@restricted <- restricted
    z@deriv <- deriv
    z@extremeScoreProbs <- c(pi.0, pi.m)
    z@rawScoresData <- rs
    z@flx.call <- z@call
    z@call <- cl
    z@nobs <-  n.total - n.0 - n.m # includes only those observations passed on to stepFlexmix
    z@identified.items <- status
  }
  else {
    ## logLik: include term for extreme scorers if estimated
    if (pi.0 != 0) z@logLiks <- z@logLiks + n.0*log(pi.0)
    if (pi.m != 0) z@logLiks <- z@logLiks + n.m*log(pi.m)
    z@models <- lapply(z@models, function(model){
      ## df: include df for extreme score probs if estimated
      model@df <- model@df + (pi.0 != 0) + (pi.m != 0)
      ## df: include df for score model if restricted over components
      if(restricted) model@df <- model@df + dfScore
      ## logLik: include term for extreme scorers if estimated
      if (pi.0 != 0) model@logLik <- model@logLik + n.0*log(pi.0)
      if (pi.m != 0) model@logLik <- model@logLik + n.m*log(pi.m)
      ## classify as raschmix and fill raschmix-specific slots
      model <- as(model, "raschmix")
      model@scores <- scores
      model@restricted <- restricted
      model@deriv <- deriv
      model@extremeScoreProbs <- c(pi.0, pi.m)
      model@rawScoresData <- rs
      model@flx.call <- model@call
      model@call <- cl
      model@call$k <- model@flx.call$k
      model@nobs <- n.total - n.0 - n.m # includes only those observations passed on to stepFlexmix
      model@identified.items <- status
      return(model)
    })
    z <- as(z, "stepRaschmix")
    z@flx.call <- z@call
    z@call <- cl
  }

  return(z)
}



## Flexmix driver for Rasch mixture model
FLXMCrasch <- function(formula = . ~ ., scores = "saturated", delta = NULL,
  nonExtremeProb = 1, ref = 1, reltol = 1e-10, deriv = "sum", hessian = FALSE, restart = TRUE, ...)
{
  scores <- match.arg(scores, c("saturated", "meanvar", "constant"))

  retval <- new("FLXMCrasch", weighted = TRUE, formula = formula, 
                name = sprintf("Rasch mixture model (%s scores)", scores))

  retval@defineComponent <- expression({
    logLik <- function(x,y,...) {
      loglikfun_rasch(item, score, y, deriv = deriv, scores = scores, nonExtremeProb = nonExtremeProb, ...)
    }
    new("FLXcomponent", df = df , logLik = logLik,
        parameters = list(item = item, score = score))
  })

  retval@fit <- function(x, y, w, component, ...){
    ## set weights to zero if computationally zero
    ## (FIXME: should this be checked by flexmix?)
    w <- ifelse(w < .Machine$double.eps^(1/1.5), 0, w)

    ## starting values
    start <- if(restart) NULL else component@parameters$item
    if (!is.null(start)){
      n <- nrow(y)
      cm <- colMeans(y[w > 0,], na.rm = TRUE)
      status <- as.character(cut(cm, c(-Inf, 1/(2 * n), 1 - 1/(2 * n), Inf), labels = c("0", "0/1", "1")))
      status[is.na(status)] <- "NA"
      status <- factor(status, levels = c("0/1", "0", "1", "NA"))
      names(status) <- colnames(y) ## FIXME: not really necessary --> remove?
      ident <- status == "0/1"
      start <- start[ident]
      start <- start[-1]
      start[!is.finite(start)] <- 0
      ## FIXME: better starting values? 0 only for NA, and very small/large value for -/+ Inf?
    }
    
    ## estimate parameters
    rasch.model <- raschmodel(y, weights = w, reltol = reltol,
      deriv = deriv, hessian = hessian, start = start, ...)

    ## # number of items and observations
    ## m <- ncol(y)
    ## #nk <- sum(w) ## obsolete

    ## ## raw scores and associated number of observations
    ## ## (some of which may be zero)
    ## rs <- as.integer(rowSums(y))
    ## nrk <- tapply(w, factor(rs, levels = 1:(m-1L)), sum)
    ## nrk[is.na(nrk)] <- 0

    ## ## raw score probabilities via logit model
    ## switch(scores,    
    ##   "saturated" = {
    ##     #fix <- min(which(nrk > 0))
    ##     xaux <- diag(m - 1L)
    ##     gamma <- log(nrk) - log(nrk[ref])
    ##   },
    ##   "meanvar" = {
    ##     rr <- as.numeric(names(nrk))
    ##     xaux <- cbind(rr / m, 4 * rr * (m - rr) / m^2)
    ##     start <- numeric(2)
    ##     clm <- function(par) {
    ##       eta <- drop(xaux %*% par)
    ##       psi <- exp(eta) / sum(exp(eta))
    ##       -sum(nrk * log(psi))
    ##     }
    ##     gamma <- optim(start, clm)$par
    ##     ref <- NULL
    ##   },      
    ##   "constant" = {
    ##     xaux <- matrix(1, nrow = m - 1L, ncol = 1L)
    ##     gamma <- 0
    ##     ref <- 1
    ##   }
    ## )

    ## ## probabilities (colSums instead of %*% to get na.rm = TRUE)
    ## eta <- colSums(as.vector(gamma) * t(xaux), na.rm = TRUE)
    ## psi <- exp(eta) / sum(exp(eta))

    ## ## drop fixed parameters
    ## if(!is.null(ref)) gamma <- gamma[-ref]
        
    ## ## any aliased parameters?
    ## alias <- if(length(gamma) > 0L) !is.finite(gamma) | is.na(gamma) else logical(0L)
    
    ## ## check for remaining problems
    ## if(any(is.na(psi))) stop("some parameters in score model not identified")
    ## ## FIXME: should this stay here?

    ## # conditional MLE -> r = 0 and r = "number of items" are excluded 
    ## #psi <- psi * nonExtremeProb ## FIXME: Do I need this here at all?

    ## ## score parameters
    ## score <- gamma#[!alias]

    ## score parameters
    if (is.null(delta)){
      scoreMod <- scoreModel(y = y, w = w, ref = ref, scores = scores)
      score <- scoreMod$score
    } else {
      score <- delta
    }
    
    ## item parameters
    cf <- c(0, rasch.model$coefficients)
    item <- structure(rep(NA, length(rasch.model$items)), .Names = names(rasch.model$items))
    item[rasch.model$items == "0/1"] <- cf
    item[rasch.model$items == "0"] <- Inf
    item[rasch.model$items == "1"] <- -Inf
    
    ## df
    df <- rasch.model$df
    if(is.null(delta)) df <- df + scoreMod$df
    
    ## collect elements necessary for loglikfun    
    #para <- list(item = worth(rasch.model), score = score, df = df)
    para <- list(item = item, score = score, df = df)

    with(para, eval(retval@defineComponent))
  }
  retval@dist <- "Rasch"
  retval
}


## compute individual contributions to log-likelihood function
loglikfun_rasch <- function(item, score, y, deriv, scores, nonExtremeProb, ...) {
  
  ## in case of non-identified parameters return log(0)  
  ##if(any(object$items != "0/1")) return(rep(-Inf, NROW(y)))
  n <- nrow(y)
  #cm <- colMeans(y, na.rm = TRUE)
  #status <- as.character(cut(cm, c(-Inf, 1/(2 * n), 1 - 1/(2 * n), Inf), labels = c("0", "0/1", "1")))
  #status[is.na(status)] <- "NA"
  #status <- factor(status, levels = c("0/1", "0", "1", "NA"))
  #ident <- status == "0/1"
  #names(status) <- colnames(y)
  #if(any(status != "0/1")) return(rep(-Inf, NROW(y)))

  ## keep orginal parameters and data
  item.o <- item
  y.o <- y
  
  ## any aliased parameters?
  alias.s <- if(length(score) > 0L) !is.finite(score) | is.na(score) else logical(0L)
  alias.i <- if(length(item) > 0L) !is.finite(item) | is.na(item) else logical(0L)
  
  ## remove aliased parameters and corresponding items (from the data)
  score <- score[!alias.s]
  item <- item[!alias.i]
  if (!isTRUE(all.equal(item[1], 0, check.attributes = FALSE))) stop("Aliased item parameter.")
  y <- y[,!alias.i]

  
  ## missing values?
  any_y_na <- any(is.na(y))
 
  ## in case of NAs
  if(any_y_na) {
    ## compute all NA patterns
    na_patterns <- factor(apply(is.na(y), 1,
                                function(z) paste(which(z), collapse = "\r")))
      
    ## replace NAs
    y[is.na(y)] <- 0
  }

  ## estimated parameters (including contrast)
  par <- item
  
  ## associated elementary symmetric functions (of order 0)
  esf <- if(any_y_na) {
    lgamma <- rep(0, NROW(y))
    for(i in levels(na_patterns)) {
      wi1 <- which(na_patterns == i)
      wi_i <- as.integer(strsplit(i, "\r")[[1]])
      par_i <- if (length(wi_i) < 1) par else par[-wi_i]
      gamma_i <- elementary_symmetric_functions(par_i,
        order = 1 - (deriv == "numeric"), diff = deriv == "diff")[[1]][-1]
      ## due to the elimination of unidentified items extreme scores are possible
      ## their loglik contributions will be set to -Inf
      lgamma[wi1] <- c(0,log(gamma_i),0)[rowSums(y[wi1, , drop = FALSE]) + 1]
    }
    lgamma
  } else {
    gamma <- elementary_symmetric_functions(par,
               order = 1 - (deriv == "numeric"), diff = deriv == "diff")[[1]][-1]
    ## due to the elimination of unidentified items extreme scores are possible
    ## their loglik contributions will be set to -Inf
    c(0,log(gamma),0)[rowSums(y)+1]
  }

  
  ## ## estimated parameters (including contrast)
  ## #par <- c(0, item) ## README: not necessary if all parameters are passed on to loglikfun_rasch
  ## par <- item
   
  ## log-likelihood contribution
  ## = parameter sum on log scale - log(sum of all permutations)
  rval <- -drop(y %*% par) - esf

  ## for the original model by Rost (1990) 
  ## log(probability of row sum r) added to log-likelihood
  ## number of items and observations
  m <- sum(colSums(y) > 0 & colSums(y) < nrow(y))
  rs <- which(tabulate(rowSums(y), ncol(y)-1L) > 0) # FIXME: use ncol(y) - 1L or m?
  nscores <- length(rs)
  switch(scores,
    "saturated" = {
      #xaux <- diag(nscores)
      xaux <- diag(length(score) + 1)
      extra <- 0
    },
    "meanvar" = {
      rr <- 1:(m-1L)
      xaux <- cbind(rr / m, 4 * rr * (m - rr) / m^2)
      extra <- NULL
    },      
    "constant" = {
      xaux <- matrix(1, nrow = m - 1L, ncol = 1L)
      extra <- 0
    }
    )

  ## calculate score probabilities
  eta <- colSums(c(extra, score) * t(xaux), na.rm = TRUE)
  psi <- exp(eta) / sum(exp(eta))
    
  ## conditional MLE -> r = 0 and r = "number of items" are excluded
  psi <- psi * nonExtremeProb

  ## check for remaining problems
  if(any(is.na(psi))) stop("some parameters in score model not identified")

  ## for saturated model:
  if (scores == "saturated"){
    ## return probability 0 for scores not present in the data
    psi.full <- numeric(m - 1L)
    psi.full[c(rs[1], as.numeric(names(score)))] <- psi
    psi <- psi.full
  }

  ## for constant model: no psi's.
  if (scores == "constant") psi <- NULL
  
  if (!is.null(psi)) rval <- c(0,log(psi),0)[rowSums(y) + 1] + rval
  ## Alternative: employ empirical distribution of the raw scores
  #if (scores == "constant") psi <- tabulate(rowSums(y), 1:(m-1))
  #rval <- log(psi)[rowSums(y)] - drop(y %*% par) - esf
  ## dann auch einfach rval in einem Schritt ausrechnen.

  ## mismatch?
  if(sum(alias.i) > 0){
    y.a <- y.o[,alias.i, drop = FALSE]
    item.a <- item.o[alias.i]
    mismatch <- matrix(NA, nrow = nrow(y.o), ncol = sum(alias.i))
    for (i in 1:sum(alias.i)){
      if(!is.na(item.a[i])){
        mismatch[,i] <- if(item.a[i] < 0) y.a[,i] == 0 else y.a[,i] == 1      
      } else {
        mismatch[,i] <- FALSE
      }
    }
    mismatch <- apply(mismatch, 1, any)
    ## set logLik contribution to smallest (machine) value possible
    ## for observation with mismatches
    rval[mismatch] <- .Machine$double.xmin
  }

  ## observations which are extreme scorers (in this iteration
  ## of the EM algorithm) are set to -Inf
  rval[rowSums(y) == 0 | rowSums(y) == ncol(y)] <- .Machine$double.xmin
  
  ## replace -Inf with smallest (machine) value possible
  ## to keep flexmix from stopping
  rval[!is.finite(rval)] <- .Machine$double.xmin

  ## return result
  return(rval)
}

mrm_clusters <- function(y, k, fit = NULL)
{
  ## data
  y <- y[, colSums(y) > 0 & colSums(y) < nrow(y), drop = FALSE]
  m <- ncol(y)
  stopifnot(m > 0)
  y <- na.omit(y)
  y <- y[rowSums(y) > 0 & rowSums(y) < m, ]

  ## raw score frequencies
  rs <- rowSums(y)
  
  ## fit mrm mixture model
  if(is.null(fit)) fit <- mRm::mrm(y, k)

  ## extract list of coefficients in RaschModel.fit parametrization
  beta <- lapply(1:k, function(i) {
    b <- fit$beta[,i]
    -b[-1L] + b[1L]
  })
    
  ## extract raw score probabilities
  psi <- lapply(1:k, function(i) fit$pi.r.c[,i])

  ## elementary symmetric functions
  lgamma <- lapply(beta, function(x) log(psychotools::elementary_symmetric_functions(c(0, x), order = 0)[[1L]][-1L]))

  ## prior probabilties
  prior <- drop(fit$class.size)
  
  ll <- sapply(1:k, function(i) -drop(y %*% c(0, beta[[i]])) - lgamma[[i]][rs] + log(psi[[i]][rs]))
  dens <- t(prior * t(exp(ll)))
  dens / rowSums(dens)
}

scoreModel <- function(y, w = NULL, ref = 1, scores = "saturated"){
  if (is.null(w)) w <- rep(1, length.out = nrow(y))
  
  # number of items and observations
  m <- ncol(y)
  
  ## raw scores and associated number of observations
  ## (some of which may be zero)
  rs <- as.integer(rowSums(y, na.rm = TRUE))
  nrk <- tapply(w, factor(rs, levels = 1:(m-1L)), sum)
  nrk[is.na(nrk)] <- 0
  ref <- min(which(nrk > 0))
  
  ## raw score probabilities via logit model
  switch(scores,    
    "saturated" = {
      #fix <- min(which(nrk > 0))
      xaux <- diag(m - 1L)
      delta <- log(nrk) - log(nrk[ref])
    },
    "meanvar" = {
      rr <- as.numeric(names(nrk))
      xaux <- cbind(rr / m, 4 * rr * (m - rr) / m^2)
      start <- numeric(2)
      clm <- function(par) {
        eta <- drop(xaux %*% par)
        psi <- exp(eta) / sum(exp(eta))
        -sum(nrk * log(psi))
      }
      delta <- optim(start, clm)$par
      ref <- NULL
    },      
    "constant" = {
      xaux <- matrix(1, nrow = m - 1L, ncol = 1L)
      delta <- 0
      ref <- 1
    }
  )
  
  ## probabilities (colSums instead of %*% to get na.rm = TRUE)
  eta <- colSums(as.vector(delta) * t(xaux), na.rm = TRUE)
  psi <- exp(eta) / sum(exp(eta))
  
  ## drop fixed parameters
  if(!is.null(ref)) delta <- delta[-ref]
  
  ## any aliased parameters?
  alias <- if(length(delta) > 0L) !is.finite(delta) | is.na(delta) else logical(0L)
  
  ## check for remaining problems
  if(any(is.na(psi))) stop("some parameters in score model not identified")
  ## FIXME: should this stay here?
  
  # conditional MLE -> r = 0 and r = "number of items" are excluded 
  #psi <- psi * nonExtremeProb ## FIXME: Do I need this here at all?

  ## score parameters
  score <- delta#[!alias]
  df <- length(delta[!alias]) ## sum(!alias)

  return(list(score = score, df = df))
}
