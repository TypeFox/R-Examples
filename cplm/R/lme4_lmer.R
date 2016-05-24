# lmer, glmer and nlmer plus methods and utilities

### Utilities for parsing the mixed model formula

findbars <- function(term)
### Return the pairs of expressions that separated by vertical bars
{
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
    if (!is.call(term)) stop("term must be of class call")
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
}

nobars <- function(term)
### Return the formula omitting the pairs of expressions that are
### separated by vertical bars
{
    if (!('|' %in% all.names(term))) return(term)
    if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
    if (length(term) == 2) {
	nb <- nobars(term[[2]])
	if (is.null(nb)) return(NULL)
	term[[2]] <- nb
	return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

subbars <- function(term)
### Substitute the '+' function for the '|' function
{
    if (is.name(term) || !is.language(term)) return(term)
    if (length(term) == 2) {
	term[[2]] <- subbars(term[[2]])
	return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name('|'))
	term[[1]] <- as.name('+')
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
}

subnms <- function(term, nlist)
### Substitute any names from nlist in term with 1
{
    if (!is.language(term)) return(term)
    if (is.name(term)) {
        if (any(unlist(lapply(nlist, get("=="), term)))) return(1)
        return(term)
    }
    stopifnot(length(term) >= 2)
    for (j in 2:length(term)) term[[j]] <- subnms(term[[j]], nlist)
    term
}

slashTerms <- function(x)
### Return the list of '/'-separated terms in an expression that
### contains slashes
{
    if (!("/" %in% all.names(x))) return(x)
    if (x[[1]] != as.name("/"))
        stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}

makeInteraction <- function(x)
### from a list of length 2 return recursive interaction terms
{
    if (length(x) < 2) return(x)
    trm1 <- makeInteraction(x[[1]])
    trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
    list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
}

expandSlash <- function(bb)
### expand any slashes in the grouping factors returned by findbars
{
    if (!is.list(bb)) return(expandSlash(list(bb)))
    ## I really do mean lapply(unlist(... - unlist returns a
    ## flattened list in this case
    unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
            return(lapply(unlist(makeInteraction(trms)),
                          function(trm) substitute(foo|bar,
                                                   list(foo = x[[2]],
                                                        bar = trm))))
        x
    }))
}

### Utilities used in lmer, glmer and nlmer

createCm <- function(A, s)
### Create the nonzero pattern for the sparse matrix Cm from A.
### ncol(A) is s * ncol(Cm).  The s groups of ncol(Cm) consecutive
### columns in A are overlaid to produce Cm.
{
    stopifnot(is(A, "dgCMatrix"))
    s <- as.integer(s)[1]
    if (s == 1L) return(A)
    if ((nc <- ncol(A)) %% s)
        stop(gettextf("ncol(A) = %d is not a multiple of s = %d",
                      nc, s))
    ncC <- as.integer(nc / s)
    TA <- as(A, "TsparseMatrix")
    as(new("dgTMatrix", Dim = c(nrow(A), ncC),
           i = TA@i, j = as.integer(TA@j %% ncC), x = TA@x),
       "CsparseMatrix")
}

### FIXME: somehow the environment of the mf formula does not have
### .globalEnv in its parent list.  example(Mmmec, package = "mlmRev")
### used to have a formula of ~ offset(log(expected)) + ... and the
### offset function was not found in eval(mf, parent.frame(2))
lmerFrames <- function(mc, formula, contrasts, vnms = character(0))
### Create the model frame, X, Y, wts, offset and terms

### mc - matched call of calling function
### formula - two-sided formula
### contrasts - contrasts argument
### vnms - names of variables to be included in the model frame
{
    mf <- mc
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]

    ## The model formula for evaluation of the model frame.  It looks
    ## like a linear model formula but includes any random effects
    ## terms and any names of parameters used in a nonlinear mixed model.
    frame.form <- subbars(formula)      # substitute `+' for `|'
    if (length(vnms) > 0)               # add the variables names for nlmer
        frame.form[[3]] <-
            substitute(foo + bar,
                       list(foo = parse(text = paste(vnms, collapse = ' + '))[[1]],
                            bar = frame.form[[3]]))

    ## The model formula for the fixed-effects terms only.
    fixed.form <- nobars(formula)       # remove any terms with `|'
    if (!inherits(fixed.form, "formula"))
      ## RHS is empty - use `y ~ 1'
      fixed.form <- as.formula(substitute(foo ~ 1, list(foo = fixed.form)))

    ## attach the correct environment
    environment(fixed.form) <- environment(frame.form) <- environment(formula)

    ## evaluate a model frame
    mf$formula <- frame.form
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fe <- mf                            # save a copy of the call
    mf <- eval(mf, parent.frame(2))

    ## evaluate the terms for the fixed-effects only (used in anova)
    fe$formula <- fixed.form
    fe <- eval(fe, parent.frame(2)) # allow model.frame to update them

    ## response vector
    Y <- model.response(mf, "any")
    ## avoid problems with 1D arrays, but keep names
    if(length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if(!is.null(nm)) names(Y) <- nm
    }
    mt <- attr(fe, "terms")

    ## Extract X checking for a null model. This check shouldn't be
    ## needed because an empty formula is changed to ~ 1 but it can't hurt.
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts) else matrix(,NROW(Y),0)
    storage.mode(X) <- "double"      # when ncol(X) == 0, X is logical
    fixef <- numeric(ncol(X))
    names(fixef) <- colnames(X)
    dimnames(X) <- NULL

    ## Extract the weights and offset.  For S4 classes we want the
    ## `not used' condition to be numeric(0) instead of NULL
    wts <- model.weights(mf); if (is.null(wts)) wts <- numeric(0)
    off <- model.offset(mf); if (is.null(off)) off <- numeric(0)

    ## check weights and offset
    if (any(wts <= 0))
        stop(gettextf("negative weights or weights of zero are not allowed"))
    if(length(off) && length(off) != NROW(Y))
        stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                      length(off), NROW(Y)))

    ## remove the terms attribute from mf
    attr(mf, "terms") <- mt
    list(Y = Y, X = X, wts = as.double(wts), off = as.double(off), mf = mf, fixef = fixef)
}

##' Is f1 nested within f2?
##'
##' Does every level of f1 occur in conjunction with exactly one level
##' of f2? The function is based on converting a triplet sparse matrix
##' to a compressed column-oriented form in which the nesting can be
##' quickly evaluated.
##'
##' @param f1 factor 1
##' @param f2 factor 2

##' @return TRUE if factor 1 is nested within factor 2

isNested <- function(f1, f2)
{
    f1 <- as.factor(f1)
    f2 <- as.factor(f2)
    stopifnot(length(f1) == length(f2))
    sm <- as(new("ngTMatrix",
                 i = as.integer(f2) - 1L,
                 j = as.integer(f1) - 1L,
                 Dim = c(length(levels(f2)),
                 length(levels(f1)))),
             "CsparseMatrix")
    all(diff(sm@p) < 2)
}


isREML <- function(x, ...) UseMethod("isREML")
isLMM  <- function(x, ...) UseMethod("isLMM")
isNLMM <- function(x, ...) UseMethod("isNLMM")
isGLMM <- function(x, ...) UseMethod("isGLMM")

##' @S3method isREML mer
isREML.mer <- function(x, ...) as.logical(x@dims["REML"])

##' @S3method isGLMM mer
isGLMM.mer <- function(x,...) {
    length(x@muEta) > 0
  ## or: is(x@resp,"glmResp")
}

##' @S3method isNLMM mer
isNLMM.mer <- function(x,...) {
  ## or: is(x@resp,"nlsResp")
  !isLMM.mer(x) & !isGLMM.mer(x)
}

##' @S3method isLMM mer
isLMM.mer <- function(x,...) as.logical(x@dims["LMM"])
## or: is(x@resp,"lmerResp") ?


##' dimsNames and devNames are in the package's namespace rather than
##' in the function lmerFactorList because the function sparseRasch
##' needs to access them.

dimsNames <- c("nt", "n", "p", "q", "s", "np", "LMM", "REML",
               "fTyp", "lTyp", "vTyp", "nest", "useSc", "nAGQ",
               "verb", "mxit", "mxfn", "cvg")
dimsDefault <- list(s = 1L,             # identity mechanistic model
                    mxit= 300L,         # maximum number of iterations
                    mxfn= 900L, # maximum number of function evaluations
                    verb= 0L,           # no verbose output
                    np= 0L,             # number of parameters in ST
                    LMM= 0L,            # not a linear mixed model
                    REML= 0L,         # glmer and nlmer don't use REML
                    fTyp= 2L,           # default family is "gaussian"
                    lTyp= 5L,           # default link is "identity"
                    vTyp= 1L, # default variance function is "constant"
                    useSc= 1L, # default is to use the scale parameter
                    nAGQ= 1L,                  # default is Laplace
                    cvg = 0L)                  # no optimization yet attempted

devNames <- c("ML", "REML", "ldL2", "ldRX2", "sigmaML",
              "sigmaREML", "pwrss", "disc", "usqr", "wrss",
              "dev", "llik", "NULLdev")


##' Create model matrices from r.e. terms.
##'
##' Create the list of model matrices from the random-effects terms in
##' the formula and the model frame.
##'
##' @param formula model formula
##' @param fr: list with '$mf': model frame; '$X': .. matrix
##' @param rmInt logical scalar - should the `(Intercept)` column
##'        be removed before creating Zt
##' @param drop logical scalar indicating if elements with numeric
##'        value 0 should be dropped from the sparse model matrices
##'
##' @return a list with components named \code{"trms"}, \code{"fl"}
##'        and \code{"dims"}
lmerFactorList <- function(formula, fr, rmInt, drop)
{
    mf <- fr$mf
    ## record dimensions and algorithm settings

    ## create factor list for the random effects
    bars <- expandSlash(findbars(formula[[3]]))
    if (!length(bars)) stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
    fl <- lapply(bars,
                 function(x)
             {
                 ff <- eval(substitute(as.factor(fac)[,drop = TRUE],
                                       list(fac = x[[3]])), mf)
                 im <- as(ff, "sparseMatrix") # transpose of indicators
		 ## Could well be that we should rather check earlier .. :
		 if(!isTRUE(validObject(im, test=TRUE)))
		     stop("invalid conditioning factor in random effect: ", format(x[[3]]))

                 mm <- model.matrix(eval(substitute(~ expr, # model matrix
                                                    list(expr = x[[2]]))),
                                    mf)
                 if (rmInt) {
                     if (is.na(icol <- match("(Intercept)", colnames(mm)))) break
                     if (ncol(mm) < 2)
                         stop("lhs of a random-effects term cannot be an intercept only")
                     mm <- mm[ , -icol , drop = FALSE]
                 }
                 ans <- list(f = ff,
                             A = do.call(rBind,
                             lapply(seq_len(ncol(mm)), function(j) im)),
                             Zt = do.call(rBind,
                             lapply(seq_len(ncol(mm)),
                                    function(j) {im@x <- mm[,j]; im})),
                             ST = matrix(0, ncol(mm), ncol(mm),
                             dimnames = list(colnames(mm), colnames(mm))))
                 if (drop) {
                     ## This is only used for nlmer models.
                     ## Need to do something more complicated for A
                     ## here.  Essentially you need to create a copy
                     ## of im for each column of mm, im@x <- mm[,j],
                     ## create the appropriate number of copies,
                     ## prepend matrices of zeros, then rBind and drop0.
                     ans$A@x <- rep(0, length(ans$A@x))
                     ans$Zt <- drop0(ans$Zt)
                 }
                 ans
             })
    dd <-
        VecFromNames(dimsNames, "integer",
                     c(list(n = nrow(mf), p = ncol(fr$X), nt = length(fl),
                            q = sum(sapply(fl, function(el) nrow(el$Zt)))),
                       dimsDefault))
    ## order terms by decreasing number of levels in the factor but don't
    ## change the order if this is already true
    nlev <- sapply(fl, function(el) length(levels(el$f)))
    ## determine the number of random effects at this point
    if (any(diff(nlev)) > 0) fl <- fl[rev(order(nlev))]
    ## separate the terms from the factor list
    trms <- lapply(fl, "[", -1)
    names(trms) <- NULL
    fl <- lapply(fl, "[[", "f")
    attr(fl, "assign") <- seq_along(fl)
    ## check for repeated factors
    fnms <- names(fl)
    if (length(fnms) > length(ufn <- unique(fnms))) {
        ## check that the lengths of the number of levels coincide
        fl <- fl[match(ufn, fnms)]
        attr(fl, "assign") <- match(fnms, ufn)
    }
    names(fl) <- ufn
    ## check for nesting of factors
    dd["nest"] <- all(sapply(seq_along(fl)[-1],
                             function(i) isNested(fl[[i-1]], fl[[i]])))
    list(trms = trms, fl = fl, dims = dd)
}

checkSTform <- function(ST, STnew)
### Check that the 'STnew' argument matches the form of ST.
{
    stopifnot(is.list(STnew), length(STnew) == length(ST),
              all.equal(names(ST), names(STnew)))
    lapply(seq_along(STnew), function (i)
           stopifnot(class(STnew[[i]]) == class(ST[[i]]),
                     all.equal(dim(STnew[[i]]), dim(ST[[i]]))))
    all(unlist(lapply(STnew, function(m) all(diag(m) > 0))))
}

lmerControl <- function(msVerbose = getOption("verbose"),
                        maxIter = 300L, maxFN = 900L)
### Control parameters for lmer, glmer and nlmer
{
    stopifnot(maxIter >= 0, maxFN >= 0)
    list(
         maxIter = as.integer(maxIter),
         maxFN = as.integer(maxFN),
	 msVerbose = as.integer(msVerbose))# "integer" on purpose
}

##' Generate a named vector of the given mode.
##' NB: If \code{defaults} contains more than one entry of a given name,
##' the *last* one wins
VecFromNames <- function(nms, mode = "numeric", defaults = list())
{
    ans <- vector(mode = mode, length = length(nms))
    names(ans) <- nms
    ans[] <- NA
    if ((nd <- length(defaults <- as.list(defaults))) > 0) {
        if (length(dnms <- names(defaults)) < nd)
            stop("defaults must be a named list")
        stopifnot(all(dnms %in% nms))
        ans[dnms] <- as(unlist(defaults), mode)
    }
    ans
}

mkZt <- function(FL, start, s = 1L)
### Create the standard versions of flist, Zt, Gp, ST, A, Cm,
### Cx, and L. Update dd.
{
    dd <- FL$dims
    fl <- FL$fl
    asgn <- attr(fl, "assign")
    trms <- FL$trms
    ST <- lapply(trms, `[[`, "ST")
    Ztl <- lapply(trms, `[[`, "Zt")
    Zt <- do.call(rBind, Ztl)
    Zt@Dimnames <- vector("list", 2)
    Gp <- c(0L, cumsum(vapply(Ztl, nrow, 1L, USE.NAMES=FALSE)))
    .Call("mer_ST_initialize", ST, Gp, Zt)
    A <- do.call(rBind, lapply(trms, `[[`, "A"))
    rm(Ztl, FL)                         # because they could be large
    nc <- sapply(ST, ncol)         # of columns in els of ST
    Cm <- createCm(A, s)
    L <- .Call("mer_create_L", Cm)
    if (s < 2) Cm <- new("dgCMatrix")
    if (!is.null(start) && checkSTform(ST, start)) ST <- start

    nvc <- sapply(nc, function (qi) (qi * (qi + 1))/2) # no. of var. comp.
### FIXME: Check number of variance components versus number of
### levels in the factor for each term. Warn or stop as appropriate

    dd["np"] <- as.integer(sum(nvc))    # number of parameters in optimization
    dev <- VecFromNames(devNames, "numeric")
    fl <- do.call(data.frame, c(fl, check.names = FALSE))
    attr(fl, "assign") <- asgn

    list(Gp = Gp, ST = ST, A = A, Cm = Cm, L = L, Zt = Zt,
         dd = dd, dev = dev, flist = fl)
}

famNms <- c("binomial", "gaussian", "Gamma", "inverse.gaussian",
            "poisson")
linkNms <- c("logit", "probit", "cauchit", "cloglog", "identity",
	     "log", "sqrt", "1/mu^2", "inverse")
varNms <- c("constant", "mu(1-mu)", "mu", "mu^2", "mu^3")

famType <- function(family)
{
    if (!(fTyp <- match(family$family, famNms, nomatch = 0)))
        stop(gettextf("unknown GLM family: %s",
                      sQuote(family$family), domain = "R-lme4"))
    if (!(lTyp <- match(family$link, linkNms, nomatch = 0)))
        stop(gettextf("unknown link: %s",
                      sQuote(family$link), domain = "R-lme4"))
    vNam <- switch(fTyp,
                   "mu(1-mu)",          # binomial
                   "constant",          # gaussian
                   "mu^2",              # Gamma
                   "mu^3",              # inverse.gaussian
                   "mu")                # poisson
    if (!(vTyp <- match(vNam, varNms, nomatch = 0)))
        stop(gettextf("unknown GLM family: %s",
                      sQuote(family$family), domain = "R-lme4"))
    c(fTyp = fTyp, lTyp = lTyp, vTyp = vTyp)
}

convergenceMessage <- function(cvg)
### Create the convergence message
{
    msg <- switch(as.character(cvg),
                  "3" = "X-convergence (3)",
                  "4" = "relative convergence (4)",
                  "5" = "both X-convergence and relative convergence (5)",
                  "6" = "absolute function convergence (6)",

                  "7" = "singular convergence (7)",
                  "8" = "false convergence (8)",
                  "9" = "function evaluation limit reached without convergence (9)",
                  "10" = "iteration limit reached without convergence (9)",
                  "14" = "storage has been allocated (?) (14)",

                  "15" = "LIV too small (15)",
                  "16" = "LV too small (16)",
                  "63" = "fn cannot be computed at initial par (63)",
                  "65" = "gr cannot be computed at initial par (65)")
    if (is.null(msg))
        msg <- paste("See PORT documentation.  Code (", cvg, ")", sep = "")
    msg
}

#### Extractors specific to mixed-effects models

coef.mer <- function(object, ...)
{
    if (length(list(...)))
        warning(paste('arguments named "',
                      paste(names(list(...)), collapse = ", "),
                      '" ignored', sep = ''))
    fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
    ref <- ranef(object)
    ## check for variables in RE but missing from FE, fill in zeros in FE accordingly
    refnames <- unlist(lapply(ref,colnames))
    nmiss <- length(missnames <- setdiff(refnames,names(fef)))
    if (nmiss >0) {
        fillvars <- setNames(data.frame(rbind(rep(0,nmiss))),missnames)
        fef <- cbind(fillvars,fef)
    }
    val <- lapply(ref, function(x) fef[rep(1, nrow(x)),,drop = FALSE])
    for (i in seq(a = val)) {
        refi <- ref[[i]]
        row.names(val[[i]]) <- row.names(refi)
        nmsi <- colnames(refi)
        if (!all(nmsi %in% names(fef)))
            stop("unable to align random and fixed effects")
        for (nm in nmsi) val[[i]][[nm]] <- val[[i]][[nm]] + refi[,nm]
    }
    class(val) <- "coef.mer"
    val
}

setMethod("coef", signature(object = "cpglmm"), coef.mer)

setAs("cpglmm", "dtCMatrix", function(from)
### Extract the L matrix
      as(from@L, "sparseMatrix"))

setMethod("fixef", signature(object = "cpglmm"),
          function(object, ...)
### Extract the fixed effects
          object@fixef)

##' Create a list of lists from multiple parallel lists

##' @param A a list
##' @param ... other, parallel lists

##' @return a list of lists

plist <- function(A, ...)
{
    dots <- list(...)
    stopifnot(is.list(A), all(sapply(dots, is.list)),
              all(sapply(dots, length) == length(A)))
    dots <- c(list(A), dots)
    ans <- A
    for (i in seq_along(A)) ans[[i]] <- lapply(dots, "[[", i)
    ans
}

##' Extract the random effects.
##'
##' Extract the conditional modes, which for a linear mixed model are
##' also the conditional means, of the random effects, given the
##' observed responses.  These also depend on the model parameters.
##'
##' @param object an object that inherits from the \code{\linkS4class{mer}} class
##' @param postVar logical scalar - should the posterior variance be returned
##' @param drop logical scalar - drop dimensions of single extent
##' @param whichel - vector of names of factors for which to return results

##' @return a named list of arrays or vectors, aligned to the factor list

setMethod("ranef", signature(object = "cpglmm"),
          function(object, postVar = FALSE, drop = FALSE, whichel = names(wt), ...)
      {
          rr <- object@ranef
          ## nt is the number of terms, cn is the list of column names
          nt <- length(cn <- lapply(object@ST, colnames))
          lterm <- lapply(plist(reinds(object@Gp), cn),
                          function(el) {
                              cni <- el[[2]]
                              matrix(rr[ el[[1]] ], ncol = length(cni),
                                     dimnames = list(NULL, cni))
                          })
          wt <- whichterms(object)
          ans <- lapply(plist(wt, object@flist),
                        function(el) {
                            ans <- do.call(cbind, lterm[ el[[1]] ])
                            rownames(ans) <- levels(el[[2]])
                            data.frame(ans, check.names = FALSE)
                        })
          ## Process whichel
          stopifnot(is(whichel, "character"))
          whchL <- names(wt) %in% whichel
          ans <- ans[whchL]

          if (postVar) {
              pV <- .Call("mer_postVar", object, whchL)
              for (i in seq_along(ans))
                  attr(ans[[i]], "postVar") <- pV[[i]]
          }
          if (drop)
              ans <- lapply(ans, function(el)
                        {
                            if (ncol(el) > 1) return(el)
                            pv <- drop(attr(el, "postVar"))
                            el <- drop(as.matrix(el))
                            if (!is.null(pv))
                                attr(el, "postVar") <- pv
                            el
                        })
          class(ans) <- "ranef.mer"
          ans
      })

print.ranef.mer <- function(x, ...) print(unclass(x), ...)
print.coef.mer <- function(x, ...) print(unclass(x), ...)

setGeneric("sigma", function(object, ...) standardGeneric("sigma"))
setMethod("sigma", signature(object = "cpglmm"),
          function (object, ...) {
              dd <- object@dims
	      if(!dd[["useSc"]]) return(1)
	      object@deviance[[if(dd[["REML"]]) "sigmaREML" else "sigmaML"]]
          })

#### Methods for standard extractors for fitted models

setMethod("anova", signature(object = "cpglmm"),
	  function(object, ...)
      {
	  mCall <- match.call(expand.dots = TRUE)
	  dots <- list(...)
	  modp <- if (length(dots))
	      sapply(dots, is, "cpglmm") | sapply(dots, is, "lm") else logical(0)
	  if (any(modp)) {		# multiple models - form table
	      opts <- dots[!modp]
	      mods <- c(list(object), dots[modp])
	      names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)],
				    as.character)
	      mods <- mods[order(sapply(lapply(mods, logLik, REML = FALSE),
					attr, "df"))]
	      calls <- lapply(mods, slot, "call")
	      data <- lapply(calls, "[[", "data")
	      if (any(data != data[[1]]))
		  stop("all models must be fit to the same data object")
	      header <- paste("Data:", data[[1]])
	      subset <- lapply(calls, "[[", "subset")
	      if (any(subset != subset[[1]]))
		  stop("all models must use the same subset")
	      if (!is.null(subset[[1]]))
		  header <-
		      c(header, paste("Subset", deparse(subset[[1]]), sep = ": "))
	      llks <- lapply(mods, logLik, REML = FALSE)
	      Df <- sapply(llks, attr, "df")
	      llk <- unlist(llks)
	      chisq <- 2 * pmax(0, c(NA, diff(llk)))
	      dfChisq <- c(NA, diff(Df))
	      val <- data.frame(Df = Df,
				AIC = sapply(llks, AIC),
				BIC = sapply(llks, BIC),
				logLik = llk,
				"Chisq" = chisq,
				"Chi Df" = dfChisq,
				"Pr(>Chisq)" = pchisq(chisq, dfChisq, lower.tail = FALSE),
				row.names = names(mods), check.names = FALSE)
	      class(val) <- c("anova", class(val))
              attr(val, "heading") <-
                  c(header, "Models:",
                    paste(rep(names(mods), times = unlist(lapply(lapply(lapply(calls,
                                           "[[", "formula"), deparse), length))),
                         unlist(lapply(lapply(calls, "[[", "formula"), deparse)),
                         sep = ": "))
	      return(val)
	  }
	  else { ## ------ single model ---------------------
            if (length(object@muEta))
              stop("single argument anova for GLMMs not yet implemented")
            if (length(object@V))
              stop("single argument anova for NLMMs not yet implemented")

            p <- object@dims[["p"]]
            ss <- (.Call("mer_update_projection", object)[[2]])^2
            names(ss) <- names(object@fixef)
            asgn <- attr(object@X, "assign")

            terms <- terms(object)
            nmeffects <- attr(terms, "term.labels")
            if ("(Intercept)" %in% names(ss))
              nmeffects <- c("(Intercept)", nmeffects)
            ss <- unlist(lapply(split(ss, asgn), sum))
            df <- unlist(lapply(split(asgn,  asgn), length))
            ## dfr <- unlist(lapply(split(dfr, asgn), function(x) x[1]))
            ms <- ss/df
            f <- ms/(sigma(object)^2)
            ## P <- pf(f, df, dfr, lower.tail = FALSE)
            ## table <- data.frame(df, ss, ms, dfr, f, P)
            table <- data.frame(df, ss, ms, f)
	    dimnames(table) <-
	      list(nmeffects,
		   ## c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
		   c("Df", "Sum Sq", "Mean Sq", "F value"))
            if ("(Intercept)" %in% nmeffects)
              table <- table[-match("(Intercept)", nmeffects), ]
            attr(table, "heading") <- "Analysis of Variance Table"
            class(table) <- c("anova", "data.frame")
            table
	  }
      })


setMethod("fitted", signature(object = "cpglmm"),
          function(object, ...)
          napredict(attr(object@frame, "na.action"), object@mu))

setMethod("residuals", signature(object = "cpglmm"),
	  function(object, ...)
          napredict(attr(object@frame, "na.action"), object@resid))

setMethod("resid", signature(object = "cpglmm"),
	  function(object, ...)
          napredict(attr(object@frame, "na.action"), object@resid))



### Show and print methods and utilities for them

formatVC <- function(varc, digits = max(3, getOption("digits") - 2))
### "format()" the 'VarCorr' matrix of the random effects -- for show()ing
{
    sc <- unname(attr(varc, "sc"))
    recorr <- lapply(varc, attr, "correlation")
    reStdDev <- c(lapply(varc, attr, "stddev"), list(Residual = sc))
    reLens <- unlist(c(lapply(reStdDev, length)))
    nr <- sum(reLens)
    reMat <- array('', c(nr, 4),
		   list(rep.int('', nr),
			c("Groups", "Name", "Variance", "Std.Dev.")))
    reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
    reMat[,2] <- c(unlist(lapply(varc, colnames)), "")
    reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
    reMat[,4] <- format(unlist(reStdDev), digits = digits)
    if (any(reLens > 1)) {
	maxlen <- max(reLens)
	corr <-
	    do.call("rBind",
		    lapply(recorr,
			   function(x, maxlen) {
			       x <- as(x, "matrix")
			       cc <- format(round(x, 3), nsmall = 3)
			       cc[!lower.tri(cc)] <- ""
			       nr <- dim(cc)[1]
			       if (nr >= maxlen) return(cc)
			       cbind(cc, matrix("", nr, maxlen-nr))
			   }, maxlen))
	colnames(corr) <- c("Corr", rep.int("", maxlen - 1))
	cbind(reMat, rBind(corr, rep.int("", ncol(corr))))
    } else reMat
}


BlockDiagonal <- function(lst)
{
    stopifnot(is(lst, "list"))
    lst <- lapply(lapply(lst, as, Class = "generalMatrix"),
                  as, Class = "TsparseMatrix")
    isSquare <- function(x) nrow(x) == ncol(x)
    stopifnot(all(sapply(lst, isSquare)),
              all(sapply(lst, is, class2 = "dMatrix")))
    if ((nl <- length(lst)) == 1) return(lst[[1]])

    offsets <- c(0L, cumsum(sapply(lst, ncol)))
    new("dgTMatrix", Dim = rep.int(offsets[nl + 1], 2),
        i = unlist(lapply(1:nl, function(i) lst[[i]]@i + offsets[i])),
        j = unlist(lapply(1:nl, function(i) lst[[i]]@j + offsets[i])),
        x = unlist(lapply(lst, slot, "x")))
}


abbrvNms <- function(gnm, cnms)
### Abbreviate names of columns in grouping factors
### gnm - group name
### cnms - column names
{
    ans <- paste(abbreviate(gnm), abbreviate(cnms), sep = '.')
    if (length(cnms) > 1) {
	anms <- lapply(cnms, abbreviate, minlength = 3)
	nmmat <- outer(anms, anms, paste, sep = '.')
	ans <- c(ans, paste(abbreviate(gnm, minlength = 3),
			    nmmat[upper.tri(nmmat)], sep = '.'))
    }
    ans
}

ST2Omega <- function(ST)
### Temporary function to convert the ST representation of the
### relative variance-covariance matrix returned by lmer into the
### Omega representation required by lmer
{
    if (nrow(ST) == 1) return(as(1/ST^2, "dpoMatrix"))
    dd <- diag(ST)
    T <- as(ST, "dtrMatrix")
    T@diag <- "U"
    crossprod(solve(T)/dd)
}


## Utilities for the fitted mer object
slotsz <- function(obj)
    rev(sort(sapply(slotNames(obj), function(s) object.size(slot(obj, s)))))

slotApply <- function(object, f, ..., simplify = FALSE) {
   .localFun <- function(what, ...) f(slot(object, what), ...)
   sapply(slotNames(object), .localFun, ..., simplify = simplify)
}


yfrm <- function(fm)
{
    stopifnot(is(fm, "cpglmm"))
    snr <- slotApply(fm, function(x)
                 {
                     if (is(x, "matrix") ||
                         is(x, "data.frame") ||
                         is(x, "numeric")) return (NROW(x))
                     0
                 }, simplify = TRUE)
    snr <- snr[snr > 0 & !(names(snr) %in%
                           c("Gp", "dims", "deviance", "frame", "flist", "X"))]
    fr <- cbind(fm@frame, fm@flist[1:NROW(fm@frame), !(names(fm@flist) %in%
                                     names(fm@frame))])
    n <- NROW(fr)
    if (NROW(fm@X) == n)
        fr <- cbind(fr, X = fm@X, Xbeta = fm@X %*% fm@fixef,
                    Zb = crossprod(fm@Zt, fm@ranef)@x)
    do.call(cbind, c(list(fr), sapply(names(which(snr == NROW(fr))),
                                      slot, object = fm, simplify = FALSE)))
}


##' Find terms associated with grouping factor names.

##' Determine the random-effects associated with particular grouping
##' factors.

##' @param fm a fitted model object of S4 class "mer"
##' @param fnm one or more grouping factor names, as a character vector

##' @return a list of indices of terms
##' @keywords models
##' @export
##' @examples
##' fm1 <- lmer(strength ~ (1|batch) + (1|sample), Pastes)
##' whichterms(fm1)
whichterms <- function(fm, fnm = names(fm@flist))
{
    stopifnot(is(fm, "cpglmm"), is.character(fnm))
    fl <- fm@flist
    asgn <- attr(fl, "assign")
    fnms <- names(fl)
    stopifnot(all(fnm %in% fnms))
    if (is.null(names(fnm))) names(fnm) <- fnm

    lapply(fnm, function(nm) which(asgn == match(nm, fnms)))
}

##' Random-effects indices by term

##' Returns a list of indices into the ranef vector by random-effects
##' terms.

##' @param Gp the Gp slot from an mer object

##' @return a list of random-effects indices
##' @keywords models
reinds <- function(Gp)
{
    lens <- diff(Gp)
    lapply(seq_along(lens), function(i) Gp[i] + seq_len(lens[i]))
}

##' Random-effects indices associated with grouping factor names

##' Determine the random-effects indices with particular grouping
##' factors.

##' @param fm a fitted model object of S4 class "mer"
##' @param fnm one or more grouping factor names, as a character vector

##' @return a list of indices of terms
##' @keywords models
##' @export
##' @examples
##' fm1 <- lmer(strength ~ (1|batch) + (1|sample), Pastes)
##' whichreind(fm1)
whichreind <- function(fm, fnm = names(fm@flist))
    lapply(whichterms(fm, fnm),
           function (ind) unlist(reinds(fm@Gp)[ind]))

