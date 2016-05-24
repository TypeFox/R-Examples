mob <- function(formula, data, subset, na.action, weights, offset, cluster,
  fit, control = mob_control(), ...)
{
  ## check fitting function
  fitargs <- names(formals(fit))
  if(!all(c("y", "x", "start", "weights", "offset") %in% fitargs)) {
    stop("no suitable fitting function specified")
  }

  ## augment fitting function (if necessary)
  if(!all(c("estfun", "object") %in% fitargs)) {
    afit <- function(y,
      x = NULL, start = NULL, weights = NULL, offset = NULL, cluster = NULL, ...,
      estfun = FALSE, object = FALSE)
    {
      obj <- if("cluster" %in% fitargs) {
        fit(y = y, x = x, start = start, weights = weights, offset = offset, cluster = cluster, ...)
      } else {
        fit(y = y, x = x, start = start, weights = weights, offset = offset, ...)
      }
      list(
        coefficients = coef(obj),
        objfun = -as.numeric(logLik(obj)),
        estfun = if(estfun) sandwich::estfun(obj) else NULL,
        object = if(object) obj else NULL
      )
    }
  } else {
    if("cluster" %in% fitargs) {
      afit <- fit
    } else {
      afit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, cluster = NULL, ..., estfun = FALSE, object = FALSE) {
        fit(y = y, x = x, start = start, weights = weights, offset = offset, ..., estfun = estfun, object = object)
      }
    }
  }

  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## formula FIXME: y ~ . or y ~ x | .
  oformula <- as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- Formula::Formula(formula(Formula::as.Formula(formula(formula), ~ 0), rhs = 2L:1L))
    xreg <- FALSE
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1L:2L))
      warning("Formula must not have more than two RHS parts")
    }
    xreg <- TRUE
  }
  mf$formula <- formula

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, response, regressor matrix (if any), partitioning variables
  mt <- terms(formula, data = data)
  mtY <- terms(formula, data = data, rhs = if(xreg) 1L else 0L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- switch(control$ytype,
    "vector" = Formula::model.part(formula, mf, lhs = 1L)[[1L]],
    "matrix" = model.matrix(~ 0 + ., Formula::model.part(formula, mf, lhs = 1L)),
    "data.frame" = Formula::model.part(formula, mf, lhs = 1L)
  )
  X <- if(!xreg) NULL else switch(control$xtype,
    "matrix" = model.matrix(mtY, mf),
    "data.frame" = Formula::model.part(formula, mf, rhs = 1L)
  )
  if(!is.null(X) && ncol(X) < 1L) {
    X <- NULL
    xreg <- FALSE
  }
  if(xreg) {
    attr(X, "formula") <- formula(formula, rhs = 1L)
    attr(X, "terms") <- mtY
    attr(X, "offset") <- cl$offset
  }
  Z <- Formula::model.part(formula, mf, rhs = 2L)
  n <- nrow(Z)
  nyx <- length(mf) - length(Z) - as.numeric("(weights)" %in% names(mf)) - as.numeric("(offset)" %in% names(mf)) - as.numeric("(cluster)" %in% names(mf))
  varindex <- match(names(Z), names(mf))

  ## weights and offset
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1L
  if(length(weights) == 1L) weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  offset <- if(xreg) model.offset(mf) else NULL
  cluster <- mf[["(cluster)"]]

  ## process pruning options (done here because of "n")
  if(!is.null(control$prune)) {
    if(is.character(control$prune)) {
      control$prune <- tolower(control$prune)
      control$prune <- match.arg(control$prune, c("aic", "bic", "none"))
      control$prune <- switch(control$prune,
        "aic" = {
	  function(objfun, df, nobs) (2 * objfun[1L] + 2 * df[1L]) < (2 * objfun[2L] + 2 * df[2L])
	}, "bic" = {
	  function(objfun, df, nobs) (2 * objfun[1L] + log(n) * df[1L]) < (2 * objfun[2L] + log(n) * df[2L])
	}, "none" = {
	  NULL
	})      
    }
    if(!is.function(control$prune)) {
      warning("Unknown specification of 'prune'")
      control$prune <- NULL
    }
  }

  ## grow the actual tree
  nodes <- mob_partynode(Y = Y, X = X, Z = Z, weights = weights, offset = offset, cluster = cluster,
    fit = afit, control = control, varindex = varindex, ...)

  ## compute terminal node number for each observation
  fitted <- fitted_node(nodes, data = mf)
  fitted <- data.frame(
      "(fitted)" = fitted,
      ## "(response)" = Y, ## probably not really needed
      check.names = FALSE,
      row.names = rownames(mf))
  if(!identical(weights, rep.int(1L, n))) fitted[["(weights)"]] <- weights
  if(!is.null(offset)) fitted[["(offset)"]] <- offset
  if(!is.null(cluster)) fitted[["(cluster)"]] <- cluster

  ## return party object
  rval <- party(nodes, 
    data = if(control$model) mf else mf[0,],
    fitted = fitted,
    terms = mt,
    info = list(
      call = cl,
      formula = oformula,
      Formula = formula,
      terms = list(response = mtY, partitioning = mtZ),
      fit = afit,
      control = control,
      dots = list(...),
      nreg = max(0L, as.integer(xreg) * (nyx - NCOL(Y))))
  )
  class(rval) <- c("modelparty", class(rval))
  return(rval)
}

## set up partynode object
mob_partynode <- function(Y, X, Z, weights = NULL, offset = NULL, cluster = NULL,
  fit, control = mob_control(), varindex = 1L:NCOL(Z), ...)
{
  ## are there regressors?
  if(missing(X)) X <- NULL
  xreg <- !is.null(X)
  n <- nrow(Z)
  if(is.null(weights)) weights <- 1L
  if(length(weights) < n) weights <- rep(weights, length.out = n)

  ## control parameters (used repeatedly)
  minsize <- control$minsize
  if(!is.null(minsize) && !is.integer(minsize)) minsize <- as.integer(minsize)
  verbose <- control$verbose
  rnam <- c("estfun", "object")
  terminal <- lapply(rnam, function(x) x %in% control$terminal)
  inner    <- lapply(rnam, function(x) x %in% control$inner)
  names(terminal) <- names(inner) <- rnam

  ## convenience functions
  w2n <- function(w) if(control$caseweights) sum(w) else sum(w > 0)
  suby <- function(y, index) {
    if(control$ytype == "vector") y[index] else y[index, , drop = FALSE]
  }
  subx <- if(xreg) {
    function(x, index) {
      sx <- x[index, , drop = FALSE]
      attr(sx, "contrasts") <- attr(x, "contrasts")
      attr(sx, "xlevels")   <- attr(x, "xlevels")
      attr(sx, "formula")   <- attr(x, "formula")
      attr(sx, "terms")     <- attr(x, "terms")
      attr(sx, "offset")    <- attr(x, "offset")
      sx
    }
  } else {
    function(x, index) NULL
  }
  subz <- function(z, index) z[index, , drop = FALSE]
  ## from strucchange
  root.matrix <- function(X) {
    if((ncol(X) == 1L)&&(nrow(X) == 1L)) return(sqrt(X)) else {
      X.eigen <- eigen(X, symmetric = TRUE)
      if(any(X.eigen$values < 0)) stop("Matrix is not positive semidefinite")
      sqomega <- sqrt(diag(X.eigen$values))
      V <- X.eigen$vectors
      return(V %*% sqomega %*% t(V))
    }
  }

  ## core mob_grow_* functions

  ## variable selection: given model scores, conduct
  ## all M-fluctuation tests for orderins in z
  mob_grow_fluctests <- function(estfun, z, weights, obj = NULL, cluster = NULL)
  {  
    ## set up return values
    m <- NCOL(z)
    pval <- rep.int(NA_real_, m)
    stat <- rep.int(0, m)
    ifac <- rep.int(FALSE, m)
    
    ## variables to test
    mtest <- if(m <= control$mtry) 1L:m else sort(sample(1L:m, control$mtry))

    ## estimating functions (dropping zero weight observations)
    process <- as.matrix(estfun)
    ww0 <- (weights > 0)
    process <- process[ww0, , drop = FALSE]
    z <- z[ww0, , drop = FALSE]
    k <- NCOL(process)
    n <- NROW(process)

    ## scale process
    process <- process/sqrt(n)
    vcov <- control$vcov
    if(is.null(obj)) vcov <- "opg"
    if(vcov != "opg") {
      bread <- vcov(obj) * n
    }
    if(vcov != "info") {
      meat <- if(is.null(cluster)) {
        crossprod(process)
      } else {
        ## nclus <- length(unique(cluster)) ## nclus / (nclus - 1L) * 
        crossprod(as.matrix(aggregate(process, by = list(cluster), FUN = sum)[, -1L, drop = FALSE]))
      }
    }
    J12 <- root.matrix(switch(vcov,
      "opg" = chol2inv(chol(meat)),
      "info" = bread,
      "sandwich" = bread %*% meat %*% bread
    ))
    process <- t(J12 %*% t(process))  

    ## select parameters to test
    if(!is.null(control$parm)) process <- process[, control$parm, drop = FALSE]
    k <- NCOL(process)

    ## get critical values for supLM statistic
    from <- if(control$trim > 1) control$trim else ceiling(n * control$trim)
    from <- max(from, minsize)
    to <- n - from
    lambda <- ((n - from) * to)/(from * (n - to))

    beta <- mob_beta_suplm
    logp.supLM <- function(x, k, lambda)
    {
      if(k > 40L) {
        ## use Estrella (2003) asymptotic approximation
        logp_estrella2003 <- function(x, k, lambda)
  	  -lgamma(k/2) + k/2 * log(x/2) - x/2 + log(abs(log(lambda) * (1 - k/x) + 2/x))
        ## FIXME: Estrella only works well for large enough x
        ## hence require x > 1.5 * k for Estrella approximation and
        ## use an ad hoc interpolation for larger p-values
        p <- ifelse(x <= 1.5 * k, (x/(1.5 * k))^sqrt(k) * logp_estrella2003(1.5 * k, k, lambda), logp_estrella2003(x, k, lambda))
      } else {
        ## use Hansen (1997) approximation
        nb <- ncol(beta) - 1L
        if(lambda<1) tau <- lambda
        else tau <- 1/(1 + sqrt(lambda))
        beta <- beta[(((k - 1) * 25 + 1):(k * 25)),]
        dummy <- beta[,(1L:nb)]%*%x^(0:(nb-1))
        dummy <- dummy*(dummy>0)
        pp <- pchisq(dummy, beta[,(nb+1)], lower.tail = FALSE, log.p = TRUE)
        if(tau == 0.5)
          p <- pchisq(x, k, lower.tail = FALSE, log.p = TRUE)
        else if(tau <= 0.01)
          p <- pp[25L]
        else if(tau >= 0.49)
          p <- log((exp(log(0.5 - tau) + pp[1L]) + exp(log(tau - 0.49) + pchisq(x, k, lower.tail = FALSE, log.p = TRUE))) * 100)
        else
        {
  	  taua <- (0.51 - tau) * 50
    	  tau1 <- floor(taua)
  	  p <- log(exp(log(tau1 + 1 - taua) + pp[tau1]) + exp(log(taua-tau1) + pp[tau1 + 1L]))
        }
      }
      return(as.vector(p))
    }

    ## compute statistic and p-value for each ordering
    for(i in mtest) {
      zi <- z[,i]
      if(is.factor(zi)) {
        proci <- process[order(zi), , drop = FALSE]
        ifac[i] <- TRUE
	iord <- is.ordered(zi) & (control$ordinal != "chisq")

        ## order partitioning variable
        zi <- zi[order(zi)]
        # re-apply factor() added to drop unused levels
        zi <- factor(zi, levels = unique(zi))
        # compute segment weights
        segweights <- as.vector(table(zi))/n

        # compute statistic only if at least two levels are left
        if(length(segweights) < 2L) {
          stat[i] <- 0
	  pval[i] <- NA_real_
        } else if(iord) {
          proci <- apply(proci, 2L, cumsum)
          tt <- head(cumsum(segweights), -1L)
          if(control$ordinal == "max") {
  	    stat[i] <- max(abs(proci[round(tt * n), ] / sqrt(tt * (1-tt))))
	    pval[i] <- log(as.numeric(1 - mvtnorm::pmvnorm(
	      lower = -stat[i], upper = stat[i],
	      mean = rep(0, length(tt)),
	      sigma = outer(tt, tt, function(x, y)
	        sqrt(pmin(x, y) * (1 - pmax(x, y)) / ((pmax(x, y) * (1 - pmin(x, y))))))
	      )^k))
	  } else {
	    proci <- rowSums(proci^2)
  	    stat[i] <- max(proci[round(tt * n)] / (tt * (1-tt)))
	    pval[i] <- log(strucchange::ordL2BB(segweights, nproc = k, nrep = control$nrep)$computePval(stat[i], nproc = k))
	  }
	} else {      
          stat[i] <- sum(sapply(1L:k, function(j) (tapply(proci[,j], zi, sum)^2)/segweights))
          pval[i] <- pchisq(stat[i], k*(length(levels(zi))-1), log.p = TRUE, lower.tail = FALSE)
        }
      } else {
        oi <- if(control$breakties) {
          mm <- sort(unique(zi))
	  mm <- ifelse(length(mm) > 1L, min(diff(mm))/10, 1)
  	  order(zi + runif(length(zi), min = -mm, max = +mm))
        } else {
          order(zi)
        }
        proci <- process[oi, , drop = FALSE]
        proci <- apply(proci, 2L, cumsum)
        stat[i] <- if(from < to) {
	  xx <- rowSums(proci^2)
	  xx <- xx[from:to]
	  tt <- (from:to)/n
	  max(xx/(tt * (1 - tt)))	  
	} else {
	  0
	}
        pval[i] <- if(from < to) logp.supLM(stat[i], k, lambda) else NA
      }
    }

    ## select variable with minimal p-value
    best <- which.min(pval)
    if(length(best) < 1L) best <- NA
    rval <- list(pval = exp(pval), stat = stat, best = best)
    names(rval$pval) <- names(z)
    names(rval$stat) <- names(z)
    if(!all(is.na(rval$best)))
      names(rval$best) <- names(z)[rval$best]
    return(rval)
  }

  ### split in variable zselect, either ordered (numeric or ordinal) or nominal
  mob_grow_findsplit <- function(y, x, zselect, weights, offset, cluster, ...)
  {
    ## process minsize (to minimal number of observations)
    if(minsize > 0.5 & minsize < 1) minsize <- 1 - minsize
    if(minsize < 0.5) minsize <- ceiling(w2n(weights) * minsize)
  
    if(is.numeric(zselect)) {
    ## for numerical variables
      uz <- sort(unique(zselect))
      if (length(uz) == 0L) stop("Cannot find admissible split point in partitioning variable")
      
      ## if starting values are not reused then the applyfun() is used for determining the split
      if(control$restart) {
        get_dev <- function(i) {
          zs <- zselect <= uz[i]
          if(w2n(weights[zs]) < minsize || w2n(weights[!zs]) < minsize) {
            return(Inf)
          } else {
            fit_left <- fit(y = suby(y, zs), x = subx(x, zs), start = NULL,
	      weights = weights[zs], offset = offset[zs], cluster = cluster[zs], ...)
            fit_right <- fit(y = suby(y, !zs), x = subx(x, !zs), start = NULL,
	      weights = weights[!zs], offset = offset[!zs], cluster = cluster[!zs], ...)
	    return(fit_left$objfun + fit_right$objfun)
          }
	}
        dev <- unlist(control$applyfun(1L:length(uz), get_dev))
      } else {
      ## alternatively use for() loop to go through all splits sequentially
      ## and reuse previous parameters as starting values
        dev <- vector(mode = "numeric", length = length(uz))

        start_left <- NULL
        start_right <- NULL

        for(i in 1L:length(uz)) {
          zs <- zselect <= uz[i]
  	  if(control$restart ||
	     !identical(names(start_left), names(start_right)) ||
	     !identical(length(start_left), length(start_right)))
	  {
	    start_left <- NULL
	    start_right <- NULL
 	  }
          if(w2n(weights[zs]) < minsize || w2n(weights[!zs]) < minsize) {
            dev[i] <- Inf
          } else {
            fit_left <- fit(y = suby(y, zs), x = subx(x, zs), start = start_left,
	      weights = weights[zs], offset = offset[zs], cluster = cluster[zs], ...)
            fit_right <- fit(y = suby(y, !zs), x = subx(x, !zs), start = start_right,
	      weights = weights[!zs], offset = offset[!zs], cluster = cluster[!zs], ...)
  	    start_left <- fit_left$coefficients
	    start_right <- fit_right$coefficients
	    dev[i] <- fit_left$objfun + fit_right$objfun
          }
        }
      }

      ## maybe none of the possible splits is admissible
      if(all(!is.finite(dev))) {
        split <- list(
          breaks = NULL,
	  index = NULL
        )
      } else {
        split <- list(
          breaks = if(control$numsplit == "center") {
	    as.double(mean(uz[which.min(dev) + 0L:1L]))
	  } else {
	    as.double(uz[which.min(dev)])
	  },
          index = NULL
        )
      }

    } else {

      if(!is.ordered(zselect) & control$catsplit == "multiway") {
        return(list(breaks = NULL, index = seq_along(levels(zselect))))      
      }

      ## for categorical variables      
      olevels <- levels(zselect) ## full set of original levels
      zselect <- factor(zselect) ## omit levels that do not occur in the data
      al <- mob_grow_getlevels(zselect)

      get_dev <- function(i) {
        w <- al[i,]
        zs <- zselect %in% levels(zselect)[w]
        if(w2n(weights[zs]) < minsize || w2n(weights[!zs]) < minsize) {
          return(Inf)
        } else {
	  if(nrow(al) == 1L) 1 else {
            fit_left <- fit(y = suby(y, zs), x = subx(x, zs), start = NULL,
	      weights = weights[zs], offset = offset[zs], cluster = cluster[zs], ...)
            fit_right <- fit(y = suby(y, !zs), x = subx(x, !zs), start = NULL,
	      weights = weights[!zs], offset = offset[!zs], cluster = cluster[zs], ...)
    	    fit_left$objfun + fit_right$objfun
	  }
        }
      }
      dev <- unlist(control$applyfun(1L:nrow(al), get_dev))

      if(all(!is.finite(dev))) {
        split <- list(
          breaks = NULL,
	  index = NULL
        )
      } else {
        if(is.ordered(zselect)) {
	  ## map back to set of full original levels
          split <- list(
            breaks = match(levels(zselect)[which.min(dev)], olevels),
	    index = NULL
          )
        } else {
	  ## map back to set of full original levels
	  ix <- structure(rep.int(NA_integer_, length(olevels)), .Names = olevels)
	  ix[colnames(al)] <- !al[which.min(dev),]
	  ix <- as.integer(ix) + 1L
          split <- list(
            breaks = NULL,
	    index = ix
          )
        }
      }
    }
  
    return(split)
  }

  ## grow tree by combining fluctuation tests for variable selection
  ## and split selection recursively
  mob_grow <- function(id = 1L, y, x, z, weights, offset, cluster, ...)
  {
    if(verbose) {
      if(id == 1L) cat("\n")
      cat(sprintf("-- Node %i %s\n", id, paste(rep("-", 32 - floor(log10(id)) + 1L), collapse = "")))
      cat(sprintf("Number of observations: %s\n", w2n(weights)))
      ## cat(sprintf("Depth: %i\n", depth))
    }

    ## fit model
    mod <- fit(y, x, weights = weights, offset = offset, cluster = cluster, ...,
      estfun = TRUE, object = terminal$object | control$vcov == "info")
    mod$test <- NULL
    mod$nobs <- w2n(weights)
    mod$p.value <- NULL

    ## set default for minsize if not specified
    if(is.null(minsize)) minsize <<- as.integer(ceiling(10L * length(mod$coefficients)/NCOL(y)))

    ## if too few observations or maximum depth: no split = return terminal node
    TERMINAL <- FALSE
    if(w2n(weights) < 2 * minsize) {
      if(verbose) cat(sprintf("Too few observations, stop splitting (minsize = %s)\n\n", minsize))
      TERMINAL <- TRUE
    }
    if(depth >= control$maxdepth) {
      if(verbose) cat(sprintf("Maximum depth reached, stop splitting (maxdepth = %s)\n\n", control$maxdepth))
      TERMINAL <- TRUE
    }
    if(TERMINAL) {
      return(partynode(id = id, info = mod))
    }

    ## conduct all parameter instability tests
    test <- if(is.null(mod$estfun)) NULL else try(mob_grow_fluctests(mod$estfun, z, weights, mod$object, cluster))

    if(!is.null(test) && !inherits(test, "try-error")) {
      if(control$bonferroni) {
        pval1 <- pmin(1, sum(!is.na(test$pval)) * test$pval)
        pval2 <- 1 - (1 - test$pval)^sum(!is.na(test$pval))
        test$pval <- ifelse(!is.na(test$pval) & (test$pval > 0.001), pval2, pval1)
      }

      best <- test$best
      TERMINAL <- is.na(best) || test$pval[best] > control$alpha
      mod$p.value <- as.numeric(test$pval[best])

      if (verbose) {
        cat("\nParameter instability tests:\n")
        print(rbind(statistic = test$stat, p.value = test$pval))
        cat(sprintf("\nBest splitting variable: %s", names(test$stat)[best]))
        cat(sprintf("\nPerform split? %s", ifelse(TERMINAL, "no\n\n", "yes\n")))
      }
    } else {
      if(verbose && inherits(test, "try-error")) cat("Parameter instability tests failed\n\n")
      TERMINAL <- TRUE
      test <- list(stat = NA, pval = NA)
    }
    
    ## update model information
    mod$test <- rbind("statistic" = test$stat, "p.value" = test$pval)

    if(TERMINAL) {
      return(partynode(id = id, info = mod))
    } else {
      zselect <- z[[best]]
      sp <- mob_grow_findsplit(y, x, zselect, weights, offset, cluster, ...)
    
      ## split successful?
      if(is.null(sp$breaks) & is.null(sp$index)) {
        if(verbose) cat(sprintf("No admissable split found in %s\n\n", sQuote(names(test$stat)[best])))
        return(partynode(id = id, info = mod))
      } else {
        sp <- partysplit(as.integer(best), breaks = sp$breaks, index = sp$index)
        if(verbose) cat(sprintf("Selected split: %s\n\n",
	  paste(character_split(sp, data = z)$levels, collapse = " | ")))
      }
    }
  
    ## actually split the data
    kidids <- kidids_split(sp, data = z)
    
    ## set-up all daugther nodes
    depth <<- depth + 1L
    kids <- vector(mode = "list", length = max(kidids))
    for(kidid in 1L:max(kidids)) {
      ## select obs for current next node
      nxt <- kidids == kidid

      ## get next node id
      if(kidid > 1L) {
        myid <- max(nodeids(kids[[kidid - 1L]]))
      } else {
        myid <- id
      }

      ## start recursion on this daugther node
      kids[[kidid]] <- mob_grow(id = myid + 1L,
        suby(y, nxt), subx(x, nxt), subz(z, nxt), weights[nxt], offset[nxt], cluster[nxt], ...)
    }
    depth <<- depth - 1L

    ## shift split varid from z to mf
    sp$varid <- as.integer(varindex[sp$varid])
    
    ## return nodes
    return(partynode(id = id, split = sp, kids = kids, info = mod))
  }

  mob_prune <- function(node)
  {
    ## turn node to list
    nd <- as.list(node)

    ## if no pruning selected
    if(is.null(control$prune)) return(nd)

    ## node information (IDs, kids, ...)
    id <- seq_along(nd)
    kids <- lapply(nd, "[[", "kids")
    tmnl <- sapply(kids, is.null)

    ## check nodes that only have terminal kids
    check <- sapply(id, function(i) !tmnl[i] && all(tmnl[kids[[i]]]))
    while(any(check)) {

      ## pruning node information
      pnode <- which(check)
      objfun <- sapply(nd, function(x) x$info$objfun)
      pok <- sapply(pnode, function(i) control$prune(
        objfun = c(objfun[i], sum(objfun[kids[[i]]])),
	df = c(length(nd[[1]]$info$coefficients), length(kids[[i]]) * length(nd[[1]]$info$coefficients) + as.integer(control$dfsplit)),
        nobs = c(nd[[i]]$info$nobs, n)
      ))

      ## do any nodes need pruning?
      pnode <- pnode[pok]
      if(length(pnode) < 1L) break

      ## prune nodes and relabel IDs  
      pkids <- sort(unlist(sapply(pnode, function(i) nd[[i]]$kids)))
      for(i in id) {
        nd[[i]]$kids <- if(nd[[i]]$id %in% pnode || is.null(kids[[i]])) {
          NULL
        } else {    
          nd[[i]]$kids - sapply(kids[[i]], function(x) sum(pkids < x))
        }
      }
      nd[pkids] <- NULL
      id <- seq_along(nd)
      for(i in id) nd[[i]]$id <- i
      
      ## node information
      kids <- lapply(nd, "[[", "kids")
      tmnl <- sapply(kids, is.null)
      check <- sapply(id, function(i) !tmnl[i] && all(tmnl[kids[[i]]]))
    }
   
    ## return pruned list 
    return(nd)
  }

  ## grow tree
  depth <- 1L
  nodes <- mob_grow(id = 1L, Y, X, Z, weights, offset, cluster, ...)

  ## prune tree
  if(verbose && !is.null(control$prune)) cat("-- Post-pruning ---------------------------\n")
  nodes <- mob_prune(nodes)  
  for(i in seq_along(nodes)) {
    if(is.null(nodes[[i]]$kids)) {
      nodes[[i]]$split <- NULL  
      if(!terminal$estfun) nodes[[i]]$info$estfun <- NULL
      if(!terminal$object) nodes[[i]]$info$object <- NULL      
    } else {
      if(!inner$estfun) nodes[[i]]$info$estfun <- NULL
      if(!inner$object) nodes[[i]]$info$object <- NULL      
    }
  }
  
  ## return as partynode
  as.partynode(nodes)
}

## determine all possible splits for a factor, both nominal and ordinal
mob_grow_getlevels <- function(z) {
  nl <- nlevels(z)
  if(inherits(z, "ordered")) {
    indx <- diag(nl)
    indx[lower.tri(indx)] <- 1
    indx <- indx[-nl, , drop = FALSE]
    rownames(indx) <- levels(z)[-nl]
  } else {
    mi <- 2^(nl - 1L) - 1L
    indx <- matrix(0, nrow = mi, ncol = nl)
    for (i in 1L:mi) {
      ii <- i
      for (l in 1L:nl) {
        indx[i, l] <- ii %% 2L
	ii <- ii %/% 2L   
      }
    }
    rownames(indx) <- apply(indx, 1L, function(x) paste(levels(z)[x > 0], collapse = "+"))
  }
  colnames(indx) <- as.character(levels(z))
  storage.mode(indx) <- "logical"
  indx
}

## control splitting parameters
mob_control <- function(alpha = 0.05, bonferroni = TRUE, minsize = NULL, maxdepth = Inf,
  mtry = Inf, trim = 0.1, breakties = FALSE, parm = NULL, dfsplit = TRUE, prune = NULL, restart = TRUE,
  verbose = FALSE, caseweights = TRUE, ytype = "vector", xtype = "matrix",
  terminal = "object", inner = terminal, model = TRUE,
  numsplit = "left", catsplit = "binary", vcov = "opg", ordinal = "chisq", nrep = 10000,
  minsplit = minsize, minbucket = minsize,
  applyfun = NULL, cores = NULL)
{
  ## transform defaults
  if(missing(minsize) & !missing(minsplit))  minsize <- minsplit
  if(missing(minsize) & !missing(minbucket)) minsize <- minbucket
  
  ## no mtry if infinite or non-positive
  if(is.finite(mtry)) {
    mtry <- if(mtry < 1L) Inf else as.integer(mtry)
  }
  
  ## data types for formula processing
  ytype <- match.arg(ytype, c("vector", "data.frame", "matrix"))
  xtype <- match.arg(xtype, c("data.frame", "matrix"))

  ## what to store in inner/terminal nodes
  if(!is.null(terminal)) terminal <- as.vector(sapply(terminal, match.arg, c("estfun", "object")))
  if(!is.null(inner))    inner    <- as.vector(sapply(inner,    match.arg, c("estfun", "object")))

  ## how to split and how to select splitting variables
  numsplit <- match.arg(tolower(numsplit), c("left", "center", "centre"))
  if(numsplit == "centre") numsplit <- "center"
  catsplit <- match.arg(tolower(catsplit), c("binary", "multiway"))
  vcov <- match.arg(tolower(vcov), c("opg", "info", "sandwich"))
  ordinal <- match.arg(tolower(ordinal), c("l2", "max", "chisq"))

  ## apply infrastructure for determining split points
  if(is.null(applyfun)) {
    applyfun <- if(is.null(cores)) {
      lapply
    } else {
      function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = cores)
    }
  }

  ## return list with all options
  rval <- list(alpha = alpha, bonferroni = bonferroni, minsize = minsize, maxdepth = maxdepth,
    mtry = mtry, trim = ifelse(is.null(trim), minsize, trim),
    breakties = breakties, parm = parm, dfsplit = dfsplit, prune = prune, restart = restart,
    verbose = verbose, caseweights = caseweights, ytype = ytype, xtype = xtype,
    terminal = terminal, inner = inner, model = model,
    numsplit = numsplit, catsplit = catsplit, vcov = vcov,
    ordinal = ordinal, nrep = nrep, applyfun = applyfun)
  return(rval)
}

## methods concerning call/formula/terms/etc.
## (default methods work for terms and update)

formula.modelparty <- function(x, extended = FALSE, ...)
  if(extended) x$info$Formula else x$info$formula

getCall.modelparty <- function(x, ...) x$info$call

model.frame.modelparty <- function(formula, ...)
{
  mf <- formula$data
  if(nrow(mf) > 0L) return(mf)

  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  mf <- formula$info$call
  mf <- mf[c(1L, match(c("formula", "data", "subset", "na.action"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf[names(nargs)] <- nargs
  if(is.null(env <- environment(formula$info$terms))) env <- parent.frame()
  mf$formula <- Formula::Formula(as.formula(mf$formula))
  eval(mf, env)
}

weights.modelparty <- function(object, ...) {
  fit <- object$fitted
  ww <- if(!is.null(w <- fit[["(weights)"]])) w else rep.int(1L, NROW(fit))
  structure(ww, .Names = rownames(fit))
}

## methods concerning model/parameters/loglik/etc.

coef.modelparty <- function(object, node = NULL, ...) {
  if(is.null(node)) node <- nodeids(object, terminal = TRUE)
  drop(do.call("rbind",
    nodeapply(object, ids = node, FUN = function(n) info_node(n)$coefficients)))
}

refit.modelparty <- function(object, node = NULL, drop = TRUE, ...)
{
  ## by default use all ids
  if(is.null(node)) node <- nodeids(object)
  
  ## estimated coefficients
  cf <- nodeapply(object, ids = node, FUN = function(n) info_node(n)$coefficients)
  
  ## model.frame
  mf <- model.frame(object)
  weights <- weights(object)
  offset <- model.offset(mf)
  cluster <- mf[["(cluster)"]]
  
  ## fitted ids
  fitted <- object$fitted[["(fitted)"]]

  ## response variables
  Y <- switch(object$info$control$ytype,
    "vector" = Formula::model.part(object$info$Formula, mf, lhs = 1L)[[1L]],
    "matrix" = model.matrix(~ 0 + ., Formula::model.part(object$info$Formula, mf, lhs = 1L)),
    "data.frame" = Formula::model.part(object$info$Formula, mf, lhs = 1L)
  )
  hasx <- object$info$nreg >= 1L | attr(object$info$terms$response, "intercept") > 0L
  X <- if(!hasx) NULL else switch(object$info$control$xtype,
    "matrix" = model.matrix(object$info$terms$response, mf),
    "data.frame" = Formula::model.part(object$info$Formula, mf, rhs = 1L)
  )
  if(!is.null(X)) {
    attr(X, "formula") <- formula(object$info$Formula, rhs = 1L)
    attr(X, "terms") <- object$info$terms$response
    attr(X, "offset") <- object$info$call$offset
  }

  suby <- function(y, index) {
    if(object$info$control$ytype == "vector") y[index] else y[index, , drop = FALSE]
  }
  subx <- if(hasx) {
    function(x, index) {
      sx <- x[index, , drop = FALSE]
      attr(sx, "contrasts") <- attr(x, "contrasts")
      attr(sx, "xlevels")   <- attr(x, "xlevels")
      attr(sx, "formula")   <- attr(x, "formula")
      attr(sx, "terms")     <- attr(x, "terms")
      attr(sx, "offset")    <- attr(x, "offset")
      sx
    }
  } else {
    function(x, index) NULL
  }
  
  ## fitting function
  afit <- object$info$fit
  
  ## refit
  rval <- lapply(seq_along(node), function(i) {
    ix <- fitted %in% nodeids(object, from = node[i], terminal = TRUE)
    args <- list(y = suby(Y, ix), x = subx(X, ix), start = cf[[i]],
      weights = weights[ix], offset = offset[ix], cluster = cluster[ix], object = TRUE)
    args <- c(args, object$info$dots)
    do.call("afit", args)$object
  })
  names(rval) <- node

  ## drop?
  if(drop & length(rval) == 1L) rval <- rval[[1L]]
  
  ## return
  return(rval)
}

apply_to_models <- function(object, node = NULL, FUN = NULL, drop = FALSE, ...) {
  if(is.null(node)) node <- nodeids(object, terminal = FALSE)
  if(is.null(FUN)) FUN <- function(object, ...) object  
  rval <- if("object" %in% object$info$control$terminal) {
    nodeapply(object, node, function(n) FUN(info_node(n)$object))
  } else {
    lapply(refit.modelparty(object, node, drop = FALSE), FUN)
  }
  names(rval) <- node
  if(drop & length(node) == 1L) rval <- rval[[1L]]
  return(rval)
}

logLik.modelparty <- function(object, dfsplit = NULL, ...)
{
  if(is.null(dfsplit)) dfsplit <- object$info$control$dfsplit
  dfsplit <- as.integer(dfsplit)
  ids <- nodeids(object, terminal = TRUE)
  ll <- apply_to_models(object, node = ids, FUN = logLik)
  dfsplit <- dfsplit * (length(object) - length(ll))
  structure(
    sum(as.numeric(ll)),
    df = sum(sapply(ll, function(x) attr(x, "df"))) + dfsplit,
    nobs = nobs(object),
    class = "logLik"
  )
}

nobs.modelparty <- function(object, ...) {
  sum(unlist(nodeapply(object,
    nodeids(object, terminal = TRUE),
    function(n) info_node(n)$nobs
  )))
}

deviance.modelparty <- function(object, ...)
{
  ids <- nodeids(object, terminal = TRUE)
  dev <- apply_to_models(object, node = ids, FUN = deviance)
  sum(unlist(dev))
}

summary.modelparty <- function(object, node = NULL, ...)
{
  ids <- if(is.null(node)) nodeids(object, terminal = TRUE) else node
  rval <- apply_to_models(object, node = ids, FUN = summary)
  if(length(ids) == 1L) rval[[1L]] else rval
}

sctest.modelparty <- function(object, node = NULL, ...)
{
  ids <- if(is.null(node)) nodeids(object, terminal = FALSE) else node
  rval <- nodeapply(object, ids, function(n) info_node(n)$test)
  names(rval) <- ids
  if(length(ids) == 1L) rval[[1L]] else rval  
}

print.modelparty <- function(x, node = NULL,
  FUN = NULL, digits = getOption("digits") - 4L,
  header = TRUE, footer = TRUE, title = NULL, objfun = "", ...)
{
  digits <- max(c(0, digits))
  if(objfun != "") objfun <- paste(" (", objfun, ")", sep = "")
  if(is.null(title)) title <- sprintf("Model-based recursive partitioning (%s)",
    deparse(x$info$call$fit))

  if(is.null(node)) {
    header_panel <- if(header) function(party) {      
      c(title, "", "Model formula:", deparse(party$info$formula), "", "Fitted party:", "")
    } else function(party) ""
  
    footer_panel <- if(footer) function(party) {
      n <- width(party)
      n <- format(c(length(party) - n, n))
      info <- nodeapply(x, ids = nodeids(x, terminal = TRUE),
        FUN = function(n) c(length(info_node(n)$coefficients), info_node(n)$objfun))
      k <- mean(sapply(info, "[", 1L))
      of <- format(sum(sapply(info, "[", 2L)), digits = getOption("digits"))

      c("", paste("Number of inner nodes:   ", n[1L]),
        paste("Number of terminal nodes:", n[2L]),
        paste("Number of parameters per node:", format(k, digits = getOption("digits"))),
        paste("Objective function", objfun, ": ", of, sep = ""), "")
    } else function (party) ""

    if(is.null(FUN)) {
      FUN <- function(x) c(sprintf(": n = %s", x$nobs), capture.output(print(x$coefficients)))
    }
    terminal_panel <- function(node) formatinfo_node(node,
      default = "*", prefix = NULL, FUN = FUN)

    print.party(x, terminal_panel = terminal_panel,
      header_panel = header_panel, footer_panel = footer_panel, ...)
  } else {
    node <- as.integer(node)
    info <- nodeapply(x, ids = node,
      FUN = function(n) info_node(n)[c("coefficients", "objfun", "test")])    
    for(i in seq_along(node)) {
      if(i == 1L) {
        cat(paste(title, "\n", collapse = ""))
      } else {
        cat("\n")
      }
      cat(sprintf("-- Node %i --\n", node[i]))
      cat("\nEstimated parameters:\n")
      print(info[[i]]$coefficients)
      cat(sprintf("\nObjective function:\n%s\n", format(info[[i]]$objfun)))
      cat("\nParameter instability tests:\n")
      print(info[[i]]$test)
    }
  }
  invisible(x)
}

predict.modelparty <- function(object, newdata = NULL, type = "node", ...)
{
  ## predicted node ids
  node <- predict.party(object, newdata = newdata)
  if(identical(type, "node")) return(node)

  ## obtain fitted model objects
  ids <- sort(unique(as.integer(node)))
  mod <- apply_to_models(object, node = ids)

  ## obtain predictions
  pred <- if(is.character(type)) {
    function(object, newdata = NULL, ...) predict(object, newdata = newdata, type = type, ...)
  } else {
    type
  }
  if("newdata" %in% names(formals(pred))) {    
    ix <- lapply(seq_along(ids), function(i) which(node == ids[i]))
    preds <- lapply(seq_along(ids), function(i)
      pred(mod[[i]], newdata = newdata[ix[[i]], , drop = FALSE], ...))
    nc <- NCOL(preds[[1L]])
    rval <- if(nc > 1L) {
      matrix(0, nrow = length(node), ncol = nc, dimnames = list(names(node), colnames(preds[[1L]])))
    } else {
      rep(preds[[1L]], length.out = length(node))
    }      
    for(i in seq_along(ids)) {
      if(nc > 1L) {
        rval[ix[[i]], ] <- preds[[i]]
	rownames(rval) <- names(node)
      } else {
        rval[ix[[i]]] <- preds[[i]]
	names(rval) <- names(node)
      }
    }
  } else {
    rval <- lapply(mod, pred, ...)
    if(NCOL(rval[[1L]]) > 1L) {
      rval <- do.call("rbind", rval)
      rownames(rval) <- ids
      rval <- rval[as.character(node), , drop = FALSE]
      rownames(rval) <- names(node)
      rval <- drop(rval)
    } else {
      ## provide a c() method for factors locally
      c.factor <- function(...) {
        args <- list(...)
	lev <- levels(args[[1L]])
	args[[1L]] <- unclass(args[[1L]])
	rval <- do.call("c", args)
	factor(rval, levels = 1L:length(lev), labels = lev)
      }
      rval <- do.call("c", rval)
      names(rval) <- ids
      rval <- rval[as.character(node)]
      names(rval) <- names(node)
    }
  }
  return(rval)
}

fitted.modelparty <- function(object, ...)
{
  ## fitted nodes
  node <- predict.party(object, type = "node")

  ## obtain fitted model objects
  ids <- nodeids(object, terminal = TRUE)
  fit <- apply_to_models(object, node = ids, FUN = fitted)

  nc <- NCOL(fit[[1L]])
  rval <- if(nc > 1L) {
    matrix(0, nrow = length(ids), ncol = nc, dimnames = list(names(node), colnames(fit[[1L]])))
  } else {
    rep(fit[[1L]], length.out = length(node))
  }	 
  for(i in seq_along(ids)) {
    if(nc > 1L) {
      rval[node == ids[i], ] <- fit[[i]]
      rownames(rval) <- names(node)
    } else {
      rval[node == ids[i]] <- fit[[i]]
      names(rval) <- names(node)
    }
  }

  return(rval)
}

residuals.modelparty <- function(object, ...)
{
  ## fitted nodes
  node <- predict.party(object, type = "node")

  ## obtain fitted model objects
  ids <- nodeids(object, terminal = TRUE)
  res <- apply_to_models(object, node = ids, FUN = residuals)

  nc <- NCOL(res[[1L]])
  rval <- if(nc > 1L) {
    matrix(0, nrow = length(ids), ncol = nc, dimnames = list(names(node), colnames(res[[1L]])))
  } else {
    rep(res[[1L]], length.out = length(node))
  }	 
  for(i in seq_along(ids)) {
    if(nc > 1L) {
      rval[node == ids[i], ] <- res[[i]]
      rownames(rval) <- names(node)
    } else {
      rval[node == ids[i]] <- res[[i]]
      names(rval) <- names(node)
    }
  }

  return(rval)
}

plot.modelparty <- function(x, terminal_panel = NULL, FUN = NULL, ...) {
  if(is.null(terminal_panel)) {
    if(is.null(FUN)) {
      FUN <- function(x) {
        cf <- x$coefficients
	cf <- matrix(cf, ncol = 1, dimnames = list(names(cf), ""))
        c(sprintf("n = %s", x$nobs), "Estimated parameters:",
          strwrap(capture.output(print(cf, digits = 4L))[-1L]))
      }
    }
    terminal_panel <- node_terminal(x, FUN = FUN)
  }
  plot.party(x, terminal_panel = terminal_panel, ...)
}
