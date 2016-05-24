varComprob <- function(fixed, data, random, groups, varcov, weights, subset, family = stats::gaussian('identity'), na.action, offset, control = varComprob.control(...), doFit = TRUE, normalizeTrace = FALSE, contrasts = NULL, model = TRUE, X = TRUE, Y = TRUE, K = TRUE, ...) {

is.formula=function(x)  
  if (missing(fixed) || !inherits(fixed, 'formula') || length(as.list(fixed))!=3L)	    stop('fixed needs to be a two-sided formula')
  if (!missing(random)) 
    .NotYetUsed("random", error = TRUE)
  if (!missing(weights)) 
    .NotYetUsed("weights", error = TRUE)
  if (!isTRUE(doFit)) { 
    .NotYetUsed("doFit", error = FALSE)
    doFit <- TRUE
  }
  if (missing(groups))
    stop('groups can not be missed')
  ## if (!missing(random)) {
  ##   if (!inherits(random, 'formula') || length(as.list(random))!=2L)
  ##     stop('random needs to be a right-sided formula, when non-missing')
  ##   if(!identical(environment(fixed) , environment(random)) )
  ##     warning("environment(random) does not match environment(fixed).")
  ##   fixedRandom <- update.formula(fixed, as.formula(paste0('.~.+',random[2],collapse=''))) # random[2] is a call to the RHS
  ##   random <- update.formula(random, as.formula('~.')) ## this will expand the formula (removing * shorthand)
  ## } else
    fixedRandom <- fixed
  if (!missing(varcov) && !is.matrix(varcov) && !is.list(varcov))
    stop("varcov needs to be a matrix or a list, when non-missing")
  if (!missing(varcov) && is.matrix(varcov))
    varcov <- list(varcov)
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (!all.equal(stats::gaussian('identity'), family)) {
    warning("Currently only Gaussian family with identity link is supported")
    family <- stats::gaussian(link='identity')
  }
  if (missing(data))
    data <- environment(formula)
  ret.x <- X; ret.y <- Y; ret.k <- K; ret.mod <- model
	
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)

  m <- match(c("fixed", "data", "random", "weights" ,"subset", "na.action", "offset", "groups"), names(mf), 0L)
  mf <- mf[c(1L, m)]  ## mf[1] is varComprob()
	
  if (m[5L]>0 && is.logical(subset)) {
    subset <- which(subset)
    mf$subset <- subset
  }
  ## if (m[4L]>0 && any(weights<=0)) { ## handling non-positive weights to subset argument
  ##   if (m[5L]>0)
  ##     subset <- unique(setdiff(subset, which(weights<=0)))
  ##   else
  ##     subset <- which(weights>0)
  ##   mf$subset <- subset
  ## }
  if (m[3L]>0L)
    mf['random'] <- NULL
  mf['fixed'] <- NULL
  mf$formula <- fixedRandom		
	
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mfAll <- eval(mf, parent.frame()) ## removing obs with missing data or not in subset
  if (!is.null( na.vec <- attr(mfAll, 'na.action'))) { ## missing data removed. 
    if (m[5L]> 0) { ## there exists subset argument
      subset <- sort(setdiff(subset, na.vec))
    } else
      subset <- sort(setdiff(seq_len(nrow(mfAll)+length(na.vec)), na.vec))
    mf$subset <- subset
    mfAll <- eval(mf, parent.frame())
  }
  mtAll <- attr(mfAll, 'terms')
  ## mf$data <- as.name('mfAll')
  mf$formula <- fixed
  mf.fixed <- eval(mf, parent.frame())	
  mt.fixed <- attr(mf.fixed, "terms")
##  fixedTerms <- attr(mt.fixed, 'term.labels')
##  fixedTerms <- sapply(fixedTerms, sortTerm)
##  fixedVars <- sapply(attr(mt.fixed, 'variables'), as.character)[-1L]
  Y <- model.response(mf.fixed, "numeric")    
  w <- as.vector(model.weights(mf.fixed))
  if (!is.null(w) && !is.numeric(w)) 
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf.fixed))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y)) 
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)", length(offset), NROW(Y)), domain = NA)
      Y <- Y - offset  ## CHECKME when family()$link is not identity
  }
  groups <- model.groups(mf.fixed)
  if (!is.null(groups)) {
    if (nrow(groups) != NROW(Y)) 
      stop(gettextf("number of groups is %d, should equal %d (number of observations)", nrow(groups), NROW(Y)), domain = NA)
  }
  if (is.empty.model(mt.fixed)) {
    X <- matrix(0, length(Y), 0L)
  } else {
    X <- model.matrix(mt.fixed, mf.fixed, contrasts)
  }
  betanames <- colnames(X)
  ## if (missing(random)) {
    if (missing(varcov))
      varcov <- list()
    K <- varcov
    nK <- length(K)
  gammanames <- names(K)
  ## } else {
  ##   ## preparing a vector of all random terms
  ##   mf$formula <- random
  ##   ## if(m[2L]>0L) mf$data <- cl$data
  ##   mf.random <- eval(mf, parent.frame())
  ##   mt.random <- attr(mf.random, 'terms')
  ##   rterms <- attr(mt.random, 'term.labels')
  ##   ## rterms <- unlist(strsplit(as.character(random[2]), ' *\\+ *'))
  ##   rterms <- setdiff(rterms, -1:1) ## removing all terms invoving intercept(s)
  ##   rterms <- sapply(rterms, sortTerm, priority=fixedVars)
  ##   rterms <- setdiff(rterms, fixedTerms)
  ##   rterms <- local({ ## expand all fixed-by-random interactions to ensure that later when fixed-by-random interactions are encounted, the marginal random effect does not need to be added. 
  ##   ans <- rterms
  ##   for (iz in seq_along(rterms)) {
  ##     this.form <- update.formula(random, as.formula(paste('~', rterms[iz], '+0', collapse='')))
  ##     mf$formula <-  this.form
  ##     mf.this <- eval(mf, parent.frame())
  ##     ## isFact <- vapply(mf.this, function(x) is.factor(x) || is.logical(x), NA)
  ##     isFixed <- names(mf.this)%in%fixedVars
  ##     ## idx <- (isFact & isFixed)
  ##     idx <- isFixed
  ##     if (any(idx) && !all(idx))
  ##       ans[iz] <- gsub(":", "*", ans[iz], fixed=TRUE)
  ##     }
  ##     ans
  ##   })
  ##   mf$formula <- update.formula(random, as.formula(paste0("~", paste0(rterms, collapse='+'))))
  ##   mf.random <- eval(mf, parent.frame())
  ##   mt.random <- attr(mf.random, 'terms')
  ##   rterms <- attr(mt.random, 'term.labels')
  ##   rterms <- setdiff(rterms, -1:1) ## removing all terms invoving intercept(s)
  ##   rterms <- sapply(rterms, sortTerm, priority=fixedVars)
  ##   rterms <- setdiff(rterms, fixedTerms)

  ##   ## processing each random term
  ##   nRterms <- length(rterms)
  ##   Z <- vector('list', nRterms)
  ##   j <- 1L
  ##   for (iz in seq_len(nRterms)) {
  ##     this.form <- update.formula(random, as.formula(paste('~', rterms[iz], '+0', collapse='')))
  ##     mf$formula <- this.form
  ##     mf.this <- eval(mf, parent.frame())
  ##     mt.this <- attr(mf.this, "terms")
  ##     tmp.model.matrix <- suppressWarnings(model.matrix(mt.this, mf.this, contrasts))
  ##     isFact <- vapply(mf.this, function(x) is.factor(x) || is.logical(x), NA)
  ##     isFixed <- names(mf.this)%in%fixedVars
  ##     idx <- which(isFact & isFixed)
  ##     if (length(idx)>0L) { ## random interactions involving fixed variables: previously, it has been ensured that random marginal effects have already been included
  ##       fterm <- paste0(sort(names(mf.this)[idx]), collapse=':')
  ##       if (fterm%in%fixedTerms) { ## truly fixed-by-random interaction
  ##         this.form <- update.formula(this.form, paste0("~",fterm))
  ##         mf$formula <- this.form
  ##         mf.tmp <- eval(mf, parent.frame())
  ##         mt.tmp <- attr(mf.tmp, 'terms')
  ##         oldContr <- getOption('contrasts')
  ##         options(contrasts=c('contr.sum','contr.sum'))
  ##         fixedX <- suppressWarnings(model.matrix(mt.tmp, mf.tmp, contrasts)[,-1L,drop=FALSE])
  ##         options(contrasts=oldContr)
  ##         tmpRterm <- substr(rterms[iz], nchar(fterm)+2L, nchar(rterms[iz]))
  ##         if (tmpRterm=='') {
  ##           warning("DEBUG ME: tmpRterm should not be empty"); browser()}
  ##         this.form <- update.formula(this.form, paste0("~", tmpRterm, '+0'))
  ##         mf$formula <- this.form
  ##         mf.tmp <- eval(mf, parent.frame())
  ##         mt.tmp <- attr(mf.tmp, 'terms')
  ##         rdmZ <- suppressWarnings(model.matrix(mt.tmp, mf.tmp, contrasts))
  ##         for (tmpi in seq_len(ncol(fixedX))) {
  ##           Z[[j]] <- fixedX[, tmpi] * rdmZ
  ##           names(Z)[j] <- paste(colnames(fixedX)[tmpi], tmpRterm, sep=':')
  ##           j <- j+1L
  ##         }
  ##       } else {  ## the fixed interaction is actually random
  ##         Z[[j]] <- tmp.model.matrix
  ##         names(Z)[j] <- rterms[iz]
  ##         j <- j+1L
  ##       }
  ##     } else {  ## no fixed terms involved
  ##       Z[[j]] <- tmp.model.matrix
  ##       names(Z)[j] <- rterms[iz]
  ##       j <- j+1L
  ##     }
  ##   }
  ##   nRterms <- length(Z)
    
  ##   if (missing(varcov)) {
  ##     nK <- nRterms
  ##     K <- vector('list', nK)
  ##     for (j in seq_len(nK))
  ##       K[[j]] <- tcrossprod(Z[[j]])
  ##   } else {
  ##     if (nRterms!=nK)
  ##       stop('The number of matrices in "varcov" needs to equal the nubmer of random effect terms in "random", when both are provided.')
  ##     for (j in 1:nK)
  ##       K[[j]] <- tcrossprod(Z[[j]]%*%K[[j]], Z[[j]])
  ##   }
  ##   names(K) <- names(Z)
  ## }
  
  if (isTRUE(normalizeTrace))
    K <- lapply(K, normalizeTrace)

  ## if (!missing(weights)) {
  ##   keep.weights <- which(w>0)
  ##   w.5 <- sqrt(w[keep.weights])
  ##   rw.5 <- w.5[rep(seq_along(keep.weights), each=length(keep.weights))]
  ##   Y <- Y[keep.weights]*w.5
  ##   X <- w.5*X[keep.weights, , drop=FALSE]
  ##   for (ik in seq_len(nK))
  ##     K[[ik]] <- w.5*K[[ik]][keep.weights, keep.weights, drop=FALSE]*rw.5
  ## }
  if (length(K)==0) {  ## linear model fit
    lm.fit(x=X, y=Y, ...) ## use lm.wfit when weights will be working!
  } else {
    g1 <- unique(groups[,1])
    g2 <- unique(groups[,2])
    p <- length(g1)
    n <- length(g2)
    k <- NCOL(X)
    YX <- arrange(data.frame(cbind(Y,X)), groups[,1], groups[,2])
    Y <- matrix(YX[,1], nrow=p, ncol=n, byrow=FALSE)
    X <- array(NA_real_, dim=c(p,n,k))
    for (j in 1:n)
      X[,j,] <- data.matrix(YX)[p*(j-1)+(1:p),-1]
    V <- array(NA_real_, dim=c(p,p,nK))
    for (j in 1:nK)
      V[,,j] <- K[[j]]
    if (!all(sapply(K, dim)==p))
      stop("the dimension of the components of the variance-covariance matrix is not the same as the one ")
    ##aperm	
    ansCall <-  call('varComprob.fit', Y=Y, X=X, V=V, control=control)
    ans <- eval(ansCall)
    names(ans$beta) <- betanames
    ans$fixef <- ans$beta
    names(ans$gamma) <- names(ans$eta) <- gammanames
    ans$parms <- ans$gamma
    names(ans$eta0) <- "error variance"
    ans$sigma2 <- ans$eta0
    ans$na.action <- list(attr(mfAll, "na.action"))  ## could be null
    ans$offset <- offset ## could be null
    ans$contrasts <- attr(X, "contrasts")  ## if X is empty, this could be null
  ## ans$xzlevels <- .getXlevels(mtAll, mfAll)
    ans$call <- cl
    ans$terms <- mtAll
    ans$nobs <- length(Y)
    ans$control <- control
    if (ret.mod)
      ans$model <- mfAll
    if (!ret.x)
      ans$X <- NULL
    else
      ans$X <- X
    if (!ret.y)
      ans$Y <- NULL
    else
      ans$Y <- Y
    if (!ret.k)
      ans$K <- NULL
    else
      ans$K <- K
    ans$random.labels <- if (is.null(names(K))) {
      if (nK > 0L)
        paste("varComp", seq_along(K), sep='.')
      else
        character(0L)
    } else
      names(K)
    ans$prior.weights <- NULL    
  ## if (missing(weights))
  ##    ans$prior.weights <- NULL
  ## else
  ##   ans$prior.weights <- w	
    class(ans) <- c(class(ans), "varComprob")
    return(ans)
  }
}
