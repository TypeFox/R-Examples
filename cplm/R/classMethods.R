
################################################
# classes defined in the cplm package
################################################

# virtual classes used in other class definitions
setClassUnion("NullNum", c("NULL","numeric"))
setClassUnion("NullList", c("NULL","list"))  
setClassUnion("NullFunc", c("NULL","function"))  
setClassUnion("ListFrame", c("list","data.frame"))

# import from package coda
setOldClass(c("mcmc", "mcmc.list", "summary.mcmc"))


## -------------------- lmer-related Classes --------------------------------

setOldClass("data.frame")
setOldClass("family")
setOldClass("logLik")

setClass("mer",
         representation(## original data
           env = "environment",# evaluation env for nonlinear model
           nlmodel = "call",# nonlinear model call
           frame = "data.frame",# model frame (or empty frame)
           call = "call",   # matched call
           flist = "data.frame",  # list of grouping factors
           X = "matrix",    # fixed effects model matrix
           Xst = "dgCMatrix", # sparse fixed effects model matrix
           Zt = "dgCMatrix",# sparse form of Z'
           pWt = "numeric",# prior weights,
           offset = "numeric", # length 0 -> no offset
           y = "numeric",   # response vector
           ###FIXME: Eliminate the cnames slot.  Put the names on the elements of the ST slot.
           #                        cnames = "list", # row/column names of els of ST
           Gp = "integer",  # pointers to row groups of Zt
           dims = "integer",# dimensions and indicators
           ## slots that vary during optimization
           ST = "list", # 
           V = "matrix",    # gradient matrix
           A = "dgCMatrix", # (ZTS)'
           Cm = "dgCMatrix", # AH'G^{-1}W^{1/2} when s > 0
           Cx = "numeric",  # x slot of Cm when s == 1 (full Cm not stored)
           L = "CHMfactor", # Cholesky factor of weighted P(AA' + I)P'
           deviance = "numeric", # ML and REML deviance and components
           fixef = "numeric",# fixed effects (length p)
           ranef = "numeric",# random effects (length q)
           u = "numeric",   # orthogonal random effects (q)
           eta = "numeric", # unbounded predictor
           mu = "numeric",  # fitted values at current beta and b
           muEta = "numeric",# d mu/d eta evaluated at current eta
           var = "numeric", # conditional variances of Y
           resid = "numeric",# raw residuals at current beta and b
           sqrtXWt = "matrix",# sqrt of model matrix row weights
           sqrtrWt = "numeric",# sqrt of weights used with residuals
           RZX = "matrix", # dense sol. to L RZX = ST'ZtX = AX
           RX = "matrix",  # Cholesky factor of downdated X'X
           ghx = "numeric", # zeros of Hermite polynomial
           ghw = "numeric"))


## -------------------- End lmer-related Classes --------------------------------

# class defining slots common to all derived classes 
setClass("cplm", 
  representation(
    call = "call",
    formula = "formula",
    contrasts = "NullList",
    link.power = "numeric",
    model.frame = "ListFrame",
    inits = "NullList")
)

# class of "cpglm", returned by a call to "cpglm" 
setClass("cpglm", 
  representation(
    coefficients = "numeric",
    residuals = "numeric",
    fitted.values = "numeric",
    linear.predictors = "numeric",
    y = "numeric",
    offset = "NullNum",
    prior.weights = "NullNum",
    weights = "numeric",
    df.residual = "integer",
    deviance = "numeric",
    aic = "numeric",
    control = "list",
    p = "numeric",
    phi = "numeric",
    iter = "integer",
    converged = "logical",
    na.action = "NullFunc",
    vcov = "matrix"),
    contains = "cplm"
)

# class of "cpglm", returned by a call to "cpglm" 
setClass("zcpglm", 
  representation(
    coefficients = "list",
    residuals = "numeric",
    fitted.values = "numeric", 
    y = "numeric",
    offset = "list",
    prior.weights = "numeric",
    df.residual = "integer",
    llik = "numeric",
    control = "list",
    p = "numeric",
    phi ="numeric",
    converged = "logical",
    na.action = "NullFunc",
    vcov = "matrix"),
    contains = "cplm"
)

# class "cpglmm" returned from a call of cpglmm
setClass("cpglmm", 
  representation(
    p = "numeric", 
    phi = "numeric",
    bound.p = "numeric",
    vcov = "matrix",
    smooths = "list"),
    contains = c("cplm", "mer")
)

# class "summary.cpglmm" 
setClass("summary.cpglmm",                 
  representation(           
    methTitle = "character",
    logLik= "logLik",
    ngrps = "integer",
    sigma = "numeric", # scale, non-negative number
    coefs = "matrix",
    REmat = "matrix",
    AICtab= "data.frame"),
  contains = "cpglmm"
)

# class "bcplm_input"
setClass("bcplm_input", 
  representation(
    X = "matrix", 
    y = "numeric",
    Zt = "dgCMatrix",
    ygt0 = "integer",
    offset = "numeric", 
    pWt = "numeric",
    mu = "numeric",
    eta = "numeric",
    inits = "list",
    fixef = "numeric",
    u = "numeric",
    phi = "numeric",
    p = "numeric",
    link.power = "numeric",
    pbeta.mean = "numeric",
    pbeta.var = "numeric",
    bound.phi = "numeric",
    bound.p = "numeric",    
    mh.sd = "numeric",              
    dims = "integer",
    k = "integer",
    Sigma = "list",
    cllik = "numeric",
    Xb = "numeric", 
    Zu = "numeric",
    Gp = "integer",
    ncol = "integer", 
    nlev = "integer",
    accept = "numeric")
)

# class of "bcplm"
setClass("bcplm", 
  representation(
    dims = "integer",  
    sims.list = "mcmc.list",
    summary = "summary.mcmc",
    prop.sd = "list",    
    Zt = "dgCMatrix",
    flist = "list",
    Sigma = "list"),
  contains="cplm"
)

         
################################################
# methods defined for cplm 
################################################

# extraction of slots using $
setMethod("$",
  signature(x = "cplm"),
  function(x, name) slot(x,name)
)

# names to get slot names
setMethod("names",
  signature(x = "cplm"),
  function(x) slotNames(x)
)

# extraction of slots using "[["
setMethod("[[",
  signature(x = "cplm", i = "numeric", j = "missing"),
  function (x, i, j, ...) slot(x,names(x)[i])
)

setMethod("[[",
  signature(x = "cplm", i = "character", j = "missing"),
  function (x, i, j, ...) slot(x, i)
)

setMethod("[",
  signature(x = "cplm", i = "numeric",
            j = "missing", drop = "missing"),
  function (x, i, j, ..., drop) {
    output <- lapply(i, function(y) slot(x, names(x)[y]))
      names(output) <- names(x)[i]
    return(output)
  }
)

setMethod("[",
  signature(x = "cplm",i = "character",
            j = "missing", drop = "missing"),
  function (x, i, j, ..., drop) {
    output <- lapply(1:length(i), function(y) slot(x, i[y]))
    names(output) <- i
    return(output)
  }
)

setMethod("terms",
  signature(x = "cplm"),
  function (x, ...) attr(x@model.frame, "terms")
)

setMethod("model.matrix",
  signature(object = "cplm"),
  function (object,...) 
    model.matrix(attr(object@model.frame, "terms"), 
          object@model.frame, object@contrasts)
)

setMethod("formula",
  signature(x = "cplm"),
  function (x, ...) x@formula
)         

setMethod("show", 
  signature(object = "cplm"),
  function(object) summary(object)                                                    
)

setMethod("vcov", 
  signature(object = "cplm"),
  function(object, ...) object@vcov
)

model.frame.cplm <- function (formula, ...) 
{
  formula@model.frame  
}
################################################
# methods defined for cpglm 
################################################
         
setMethod("coef",
  signature(object = "cpglm"),
  function (object, ...) object@coefficients
)

setMethod("residuals",
  signature(object = "cpglm"),
  function (object, type = c("deviance", "pearson", "working", 
    "response", "partial"), ...) {      
    type <- match.arg(type)
    y <- object@y
    r <- object@residuals
    mu <- object@fitted.values
    wts <- object@prior.weights
    family <- tweedie(var.power = object@p,link.power = object@link.power)
    switch(type, deviance = , pearson = , response = if (is.null(y)) {
        eta <- object@linear.predictors
        y <- mu + r * family$mu.eta(eta)
    })
    res <- switch(type, 
      deviance = if (object@df.residual > 0) {
        d.res <- sqrt(pmax((family$dev.resids)(y, mu, 
            wts), 0))
        ifelse(y > mu, d.res, -d.res)
        } else rep.int(0, length(mu)), 
      pearson = (y - mu) * sqrt(wts)/sqrt(family$variance(mu)), 
      working = r, 
      response = y - mu, 
      partial = r)
    na.action <- attr(object@model.frame,"na.action")
    if (!is.null(na.action)) 
        res <- naresid(na.action, res)
    #if (type == "partial") 
    #    res <- res + predict(object, type = "terms")
    res
    }
)

setMethod("resid",
  signature(object = "cpglm"),
  function (object, type = c("deviance", "pearson", "working", 
  "response", "partial"), ...) 
    return(residuals(object, type = type))
)

# generate fitted values on the original scale
setMethod("fitted",
  signature(object = "cpglm"),
  function (object, ...) object@fitted.values
)
	
setMethod("AIC",
  signature(object = "cpglm",k = "missing" ),
  function (object, ..., k) object@aic
)

setMethod("deviance",
  signature(object = "cpglm"),
  function (object, ...) object@deviance
)

setMethod("summary", signature(object = "cpglm"),
	function(object,...){
    coef.beta <- coef(object)  
    vc <- vcov(object)
    s.err <- sqrt(diag(vc))    
    err.beta <- s.err
    test.value <- coef.beta / err.beta
    dn <- c("Estimate", "Std. Error")             
    pvalue <- 2 * pt(-abs(test.value), object@df.residual)
    
    coef.table <- cbind(coef.beta, err.beta, test.value, pvalue)  
    dn2 <- c("t value", "Pr(>|t|)")
    dimnames(coef.table) <- list(names(coef.beta), c(dn, dn2))
    keep <- match(c("call", "deviance", "aic", "contrasts", "df.residual",  
        "iter","na.action"), names(object), 0L)  
    ans <- c(object[keep], list(deviance.resid = residuals(object, 
        type = "deviance"), coefficients = coef.table, 
        dispersion = object@phi, vcov = vc, p = object@p))    
    .print.cpglm.summary(ans)    
    }
)

.print.cpglm.summary<-function(x,digits = max(3, getOption("digits") - 3),
                               signif.stars = getOption("show.signif.stars"), ...){
  
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("Deviance Residuals: \n")
    if (x$df.residual > 5) {
        x$deviance.resid <- quantile(x$deviance.resid, na.rm = TRUE)
        names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q", 
            "Max")
    }
    xx <- zapsmall(x$deviance.resid, digits + 1)
    print.default(xx, digits = digits, na.print = "", print.gap = 2)
    printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, 
            na.print = "NA",...)
        
    cat("\nEstimated dispersion parameter:",  
        format(x$dispersion, digits = max(5, digits + 1))) 
    cat("\nEstimated index parameter:",  
        format(x$p, digits = max(5, digits + 1)),"\n\n") 
    cat("Residual deviance:", format(x$deviance, digits = max(5, digits + 1)), 
        " on", format(x$df.residual), " degrees of freedom\n") 
    if (nzchar(mess <- naprint(x$na.action))) 
        cat("  (", mess, ")\n", sep = "")
    cat("AIC: ", format(x$aic, digits = max(4, digits + 1)), "\n\n")
    cat("Number of Fisher Scoring iterations: ", x$iter, "\n") 
    cat("\n")
    invisible(x)
}

# simple prediction method for cpglm
setMethod("predict", signature(object = "cpglm"),
  function (object, newdata, type = c("response", "link"), 
                  na.action = na.pass, ...) {
    tt <- attr(object@model.frame, "terms")
    if (missing(newdata) || is.null(newdata)) {
        X <- model.matrix(object)
        offset <- object$offset
    }
    else {
        Terms <- delete.response(tt)
        xlevels <- .getXlevels(Terms, object@model.frame)
        m <- model.frame(Terms, newdata, na.action = na.action, xlev = xlevels)
        X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
        offset <- rep(0, nrow(X))
        if (!is.null(off.num <- attr(tt, "offset"))) 
            for (i in off.num) offset <- offset + eval(attr(tt, 
                "variables")[[i + 1]], newdata)
        if (!is.null(object$call$offset)) 
            offset <- offset + eval(object$call$offset, newdata)
    }
   beta <- object$coefficients
   na.ps <- which(is.na(beta))
   if (length(na.ps)) {
    predictor <- X[, -na.ps, drop = FALSE] %*% beta[-na.ps]
    warning("prediction from a rank-deficient fit may be misleading")
   } else {
    predictor <- X%*% beta
   }
    if (!is.null(offset)) 
        predictor <- predictor + offset
    mu <- tweedie(link.power = object@link.power)$linkinv(predictor)
    type <- match.arg(type)                                                            
    switch(type,link = predictor, response = mu)                                                            
})


        
################################################
# methods defined for cpglmm
################################################


setMethod("vcov", signature(object = "cpglmm"),
  function(object, ...){
    rr <- object$phi * chol2inv(object@RX, size = object@dims['p'])
    nms <- colnames(object@X)
    dimnames(rr) <- list(nms, nms)
    if (FALSE){  
      # compute vcov for phi and p numerically 
      cpglmm_dev <- function(x, ...){
        parm <- c(.Call("cpglmm_ST_getPars", object), 
                  object$fixef, log(x[1]), x[2])
        .Call("cpglmm_update_dev", object, parm)  
      }
      x <- c(object$phi, object$p)
      hs <- hess(x, cpglmm_dev)
      dimnames(hs) <- list(c("phi", "p"), c("phi", "p"))      
      attr(rr,"phi_p") <- solve(hs)
    }
    rr
})

setGeneric("VarCorr", function(x, ...) standardGeneric("VarCorr"))
setMethod("VarCorr", signature(x = "cpglmm"),
  function(x, ...){
    sc <- sqrt(x@phi)
    ans <- lapply(cc <- .Call("cpglmm_ST_chol", x),
             function(ch) {
                val <- crossprod(sc * ch) # variance-covariance
                stddev <- sqrt(diag(val))
                correl <- t(val / stddev)/stddev
                diag(correl) <- 1
                attr(val, "stddev") <- stddev
                attr(val, "correlation") <- correl
                val
              })
    fl <- x@flist
    names(ans) <- names(fl)[attr(fl, "assign")]
    attr(ans, "sc") <- sc
    ans
})



setMethod("logLik", signature(object="cpglmm"),
          function(object, REML = NULL, ...)
            ### Extract the log-likelihood or restricted log-likelihood
          {
            dims <- object@dims
            if (is.null(REML) || is.na(REML[1]))
              REML <- dims[["REML"]]
            val <- -object@deviance["ML"]/2
            attr(val, "nall") <- attr(val, "nobs") <- dims[["n"]]
            attr(val, "df") <-
              dims[["p"]] + dims[["np"]] + as.logical(dims[["useSc"]])
            attr(val, "REML") <-  as.logical(REML)
            class(val) <- "logLik"
            val
          })

setMethod("summary", signature(object = "cpglmm"),
  function(object, ...){
    fcoef <- fixef(object)
    vcov <- object@vcov
    dims <- object@dims
    coefs <- cbind("Estimate" = fcoef, "Std. Error" = sqrt(diag(vcov)) )
    llik <- logLik(object)
    dev <- object@deviance
    mType <- "LMM"
    mName <- "Compound Poisson linear"
    method <- paste("the", if(dims[["nAGQ"]] == 1) "Laplace" else
		  "adaptive Gaussian Hermite","approximation")
  
    AICframe <- data.frame(AIC = AIC(llik), BIC = BIC(llik),
                               logLik = as.vector(llik),
                               deviance = dev[["ML"]],
                               row.names = "")
    
    varcor <- VarCorr(object)
    REmat <- formatVC(varcor)
    if (is.na(attr(varcor, "sc")))
        REmat <- REmat[-nrow(REmat), , drop = FALSE]

    if (nrow(coefs) > 0) {
      if (!dims[["useSc"]]) {
        coefs <- coefs[, 1:2, drop = FALSE]
        stat <- coefs[,1]/coefs[,2]
        pval <- 2*pnorm(abs(stat), lower.tail = FALSE)
        coefs <- cbind(coefs, "z value" = stat, "Pr(>|z|)" = pval)
      } else {
        stat <- coefs[,1]/coefs[,2]
        ##pval <- 2*pt(abs(stat), coefs[,3], lower = FALSE)
        coefs <- cbind(coefs, "t value" = stat) #, "Pr(>|t|)" = pval)
      }
    } 
    new("summary.cpglmm", object,
            methTitle = paste(mName, "mixed model fit by", method),
            logLik = llik,
            ngrps = sapply(object@flist, function(x) length(levels(x))),
            sigma = sqrt(object@phi),
            coefs = coefs,
            REmat = REmat,
            AICtab = AICframe)
  }
)

## This is modeled a bit after  print.summary.lm :
print.cpglmm <- function(x, digits = max(3, getOption("digits") - 3),
                     correlation = FALSE, symbolic.cor = FALSE,
                     signif.stars = getOption("show.signif.stars"), ...){
  so <- summary(x)
  llik <- so@logLik
  dev <- so@deviance
  dims <- x@dims
  cat(so@methTitle, "\n")
  if (!is.null(x@call$formula))
      cat("Formula:", deparse(x@call$formula),"\n")
  if (!is.null(x@call$data))
      cat("   Data:", deparse(x@call$data), "\n")
  if (!is.null(x@call$subset))
      cat(" Subset:", deparse(x@call$subset),"\n")
  print(so@AICtab, digits = digits)

  cat("Random effects:\n")
  print(so@REmat, quote = FALSE, digits = digits, ...)

  ngrps <- so@ngrps
  cat(sprintf("Number of obs: %d, groups: ", dims[["n"]]))
  cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
  cat("\n")
   
  if (nrow(so@coefs) > 0) {
    cat("\nFixed effects:\n")
    printCoefmat(so@coefs, zap.ind = 3, #, tst.ind = 4
	     digits = digits, signif.stars = signif.stars)

    cat("\nEstimated dispersion parameter:", round(so@phi, digits=digits))
    cat("\n")
    cat("Estimated index parameter:", round(so@p, digits=digits))
    cat("\n")
     
    if(correlation) {
      corF <- so@vcov@factors$correlation
      if (!is.null(corF)) {
  	    p <- ncol(corF)
  	    if (p > 1) {
  	      rn <- rownames(so@coefs)
  	      rns <- abbreviate(rn, minlength=11)
  	      cat("\nCorrelation of Fixed Effects:\n")
  	      if (is.logical(symbolic.cor) && symbolic.cor) {
  		      corf <- as(corF, "matrix")
  		      dimnames(corf) <- list(rns,
  				       abbreviate(rn, minlength=1, strict=TRUE))
  		      print(symnum(corf))
  	      } else {
  		      corf <- matrix(format(round(corF@x, 3), nsmall = 3),
  			       ncol = p,dimnames = list(rns, abbreviate(rn, minlength=6)))
  		      corf[!lower.tri(corf)] <- ""
  		      print(corf[-1, -p, drop=FALSE], quote = FALSE)
  	      }
  	    }
      }
    }
  }
  invisible(x)
}

setMethod("print", "cpglmm", print.cpglmm)
setMethod("show", "cpglmm", 
  function(object) print.cpglmm(object)
)


# predict method for cpglmm
getZt <- function(formula, oldmf, newmf){
  bars <- expandSlash(findbars(formula[[3]]))
  names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
  fl <- lapply(bars, function(x) {
        oldlvl <- eval(substitute(levels(as.factor(fac)[, drop = TRUE]), 
            list(fac = x[[3]])), oldmf)
        ff <- eval(substitute(factor(fac,levels = oldlvl)[, drop = TRUE], 
            list(fac = x[[3]])), newmf)
        # fill columns of 0's if some levels are missing
        im <- as(ff, "sparseMatrix")
        im2 <- Matrix(0, nrow = length(oldlvl), ncol = length(ff), sparse = TRUE)
        # this is awkward as the Matrix package seems to fail
        for (i in 1:nrow(im)){
          ind <- match(rownames(im)[i], oldlvl)
          im2[as.numeric(ind), ] <- im[as.numeric(i), ]            
        }        
        if (!isTRUE(validObject(im, test = TRUE))) 
            stop("invalid conditioning factor in random effect: ", 
                format(x[[3]]))
        mm <- model.matrix(eval(substitute(~expr, list(expr = x[[2]]))), newmf)
        mm <- mm[!is.na(ff), , drop = F]
        Zt <- do.call(rBind, lapply(seq_len(ncol(mm)), 
            function(j) {
                im2@x <- mm[, j]
                im2
            }))
        ans <- list(f = oldlvl, Zt = Zt)       
        ans
    })
  nlev <- sapply(fl, function(el) length(levels(el$f)))
  if (any(diff(nlev)) > 0) 
        fl <- fl[rev(order(nlev))]        
  Zt <- do.call(rBind, lapply(fl, "[[", "Zt"))
  Zt
}         

setMethod("predict", signature(object = "cpglmm"),
    function(object,  newdata, type = c("response", "link"), 
          na.action = na.pass, ...) {
    tt <- attr(object@model.frame,"terms")
    if (missing(newdata) || is.null(newdata)) {
        mm <- X <- model.matrix(object)
        Zt <- object@Zt
        offset <- object$offset
    }
    else {
      #FIXME: should I use xlev ???
        Terms <- delete.response(tt)
        # design matrix for fixed effects
        X <- model.matrix(Terms, newdata, contrasts.arg = object@contrasts)
        # design matrix for random effects
        formula <- object@formula
        oldmf <- object@model.frame
        Zt <- getZt(formula, oldmf, newdata)        
        # get offset
        offset <- rep(0, nrow(X))
        if (!is.null(off.num <- attr(tt, "offset"))) 
            for (i in off.num) offset <- offset + eval(attr(tt, 
                "variables")[[i + 1]], newdata)
        if (!is.null(object$call$offset)) 
            offset <- offset + eval(object$call$offset, newdata)                
    }
    beta <- object@fixef
    u <- object@ranef
    predictor <- as.numeric(X %*% beta + t(Zt)%*% u)
    if (!is.null(offset)) 
        predictor <- predictor + offset
    mu <- tweedie(link.power = object@link.power)$linkinv(predictor)
    type <- match.arg(type)
    switch(type,link = predictor, response = mu)   
})

################################################
# methods defined for bcplm
################################################

# fixed effects
setMethod("fixef", signature = "bcplm",
  function(object, type = c("median", "mean"), sd = FALSE,
           quantiles = NULL, ...){
    type <- match.arg(type)
    s <- object@summary
    dm <- object@dims
    rn <- 1:unname(dm["n.beta"])
    mu.beta <- if (type == "median") as.numeric(s[[2]][rn, 3]) else 
      as.numeric(s[[1]][rn, 1])
    names(mu.beta) <- rownames(s[[1]])[rn]
    if (sd){
      sd.beta <- as.numeric(s[[1]][rn, 2])
      attr(mu.beta, "sd") <- sd.beta
    }
    if (!is.null(quantiles)){
      qt <- as.matrix(summary(object$sims.list, quantiles = quantiles)[[2]])
      attr(mu.beta, "quantiles") <- qt[rn, , drop = FALSE]
    }
    return(mu.beta)
  }
)

# variance components
setMethod("VarCorr", signature(x = "bcplm"),
  function(x, ...){
    dm <- x@dims
    if (dm["n.u"] == 0)
      stop("No random effects in 'VarCorr'!")
    
    ans <- lapply(x@Sigma, function(xx) {
      stddev <- sqrt(diag(xx))
      correl <- t(xx / stddev)/stddev
      diag(correl) <- 1
      attr(xx, "stddev") <- stddev
      attr(xx, "correlation") <- correl
      xx
    })
    fl <- x@flist
    names(ans) <- names(fl)[attr(fl, "assign")]
    attr(ans, "sc") <- sqrt(x@summary[[2]][dm["n.beta"] + 1, 3])
    ans
  }
)

setMethod("show", signature = "bcplm",  
  function(object) 
    print.bcplm(object)
)

setMethod("summary", signature = "bcplm",  
  function(object) 
    object
)

setMethod("plot", signature(x = "bcplm", y = "missing"),
  function(x, y, ...) plot(x@sims.list)
)


# print out (summarize) model results
print.bcplm <- function(x, digits = max(3, getOption("digits") - 3)){
  dims <- x@dims
  # fixed effects
  fcoef <- fixef(x, sd = TRUE, quantiles = c(0.025, 0.975))
  coefs <- cbind("Estimate" = fcoef, "Std. Error" = attr(fcoef, "sd"),
                 "Lower (2.5%)" = attr(fcoef, "quantiles")[, 1],
                 "Upper (97.5%)" = attr(fcoef, "quantiles")[, 2])
  
  # start printing
  cat("Compound Poisson linear models via MCMC\n")
  cat(dims["n.chains"], " chains, each with ", dims["n.iter"], " iterations (first ", 
      dims["n.burnin"], " discarded)", sep = "")
  if (dims["n.thin"] > 1) 
    cat(", n.thin =", dims["n.thin"])
  cat("\nn.sims =", dims["n.sims"], "iterations saved\n")
  cat("\n")
  if (!is.null(x@call$formula))
    cat("Formula:", deparse(x@call$formula),"\n")
  if (!is.null(x@call$data))
    cat("   Data:", deparse(x@call$data), "\n")
  if (!is.null(x@call$subset))
    cat(" Subset:", deparse(x@call$subset),"\n")
  
  if (dims["n.u"] > 0){
    cat("\nRandom and dynamic variance components:\n")
    varcor <- VarCorr(x)
    REmat <- formatVC(varcor)
    if (is.na(attr(varcor, "sc")))
      REmat <- REmat[-nrow(REmat), , drop = FALSE]      
    print(REmat, quote = FALSE, digits = digits)
  }
  cat(sprintf("Number of obs: %d ", x@dims["n.obs"]))
  if (dims["n.u"] > 0){
    ngrps <- sapply(x@flist, nlevels)
    cat(", groups: ")
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
  }
  cat("\n")
  
  if (nrow(coefs) > 0) {
    cat("\nFixed effects:\n")
    printCoefmat(coefs, zap.ind = 3, digits = digits)  
    cat("---")
  }
  s <- x@summary
  phi.ps <- grep("^phi$", rownames(s[[1]]))
  p.ps <- grep("^p$", rownames(s[[1]]))
  cat("\nEstimated dispersion parameter:",  
      format(s[[2]][phi.ps, 3], digits = max(5, digits + 1))) 
  cat("\nEstimated index parameter:",  
      format(s[[2]][p.ps, 3], digits = max(5, digits + 1)),"\n\n")
  out <- list(fixef = coefs, 
              VarCorr = if (dims["n.u"]) REmat else list())
  invisible(out)
}


################################################
# methods defined for zcpglm 
################################################
         
setMethod("coef",
  signature(object = "zcpglm"),
  function (object, ...) object@coefficients
)

setMethod("residuals",
  signature(object = "zcpglm"), 
  function(object, ...) object@residuals  
)

setMethod("resid",
  signature(object = "zcpglm"), 
  function(object, ...) residuals(object)
)

# generate fitted values on the original scale
setMethod("fitted",
  signature(object = "zcpglm"),
  function (object, ...)  object@fitted.values
)
	
setMethod("summary", signature(object = "zcpglm"),
	function(object, ...){ 
    nbz <- length(coef(object)$zero)
    nbt <- length(coef(object)$tweedie)
    se <- sqrt(diag(vcov(object)))
    coef <- unlist(coef(object))
    zstat <- coef / se
    pval <- 2 * pnorm(-abs(zstat))
    coef <- cbind(coef, se, zstat, pval)
    colnames(coef) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(coef) <- c(names(coef(object)$zero), names(coef(object)$tweedie))
    
    coef.table <- list()
    coef.table$zero <- coef[1:nbz, , drop = FALSE]
    coef.table$tweedie <- coef[(nbz + 1):(nbz + nbt), , drop = FALSE]
    
    keep <- match(c("llik", "contrasts", "df.residual",  
        "na.action", "vcov"), names(object), 0L)
    
    out <- list(llik = object@llik, contrasts = object@contrasts, 
                df.residual = object@df.residual, vcov = object @vcov, 
                na.action = object@na.action, coefficients = coef.table,
                call = object@call, phi = object@phi, p = object@p)
    .print.zcpglm.summary(out)  
	}
)
                    
.print.zcpglm.summary<-function(x,digits = max(3, getOption("digits") - 3),
                               signif.stars = getOption("show.signif.stars"), ...){
  
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat(paste("Zero-inflation model coefficients:\n"))
    printCoefmat(x$coefficients$zero, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", signif.legend = FALSE)
    cat(paste("\nCompound Poisson model coefficients:\n"))
    printCoefmat(x$coefficients$tweedie, digits = digits, signif.stars = signif.stars, 
            na.print = "NA")    
    cat("\nEstimated dispersion parameter:",  
        format(x$phi, digits = max(5, digits + 1))) 
    cat("\nEstimated index parameter:",  
        format(x$p, digits = max(5, digits + 1)),"\n") 
    if (nzchar(mess <- naprint(x$na.action))) 
        cat("  (", mess, ")\n", sep = "")
    invisible(x)
}


setMethod("predict", signature(object = "zcpglm"),
  function (object, newdata, type = c("response", "zero", "tweedie"), 
                    na.action = na.pass, ...) {
  call <- object$call
  ttz <- attr(object@model.frame$zero, "terms")
  ttt <- attr(object@model.frame$tweedie, "terms")
  Termsz <- delete.response(ttz)
  Termst <- delete.response(ttt)
  xlevz <- .getXlevels(Termsz, object@model.frame$zero)
  xlevt <- .getXlevels(Termst, object@model.frame$tweedie)  
  mz <- model.frame(Termsz, newdata, na.action = na.action, xlev = xlevz)
  mt <- model.frame(Termst, newdata, na.action = na.action, xlev = xlevt)
  Xz <- model.matrix(Termsz, mz, contrasts.arg = object$contrasts)
  Xt <- model.matrix(Termst, mt, contrasts.arg = object$contrasts)
  offt <- offz <- rep(0, nrow(Xz))
  if (!is.null(off.num <- attr(ttz, "offset"))) 
    for (i in off.num) 
        offz <- offz + eval(attr(ttz, "variables")[[i + 1]], newdata)
  if (!is.null(off.num <- attr(ttt, "offset"))) 
    for (i in off.num) 
      offt <- offt + eval(attr(ttt, "variables")[[i + 1]], newdata)
  if (!is.null(object$call$offset)) {
      off <- eval(object$call$offset, newdata)
      offz <- offz + off
      offt <- offt + off 
  }
  
  link.power <- make.link.power(object$link.power)
  tw <- tweedie(link.power = link.power)                 
  logit <- binomial()
  
  betaz <- object$coefficients$zero
  betat <- object$coefficients$tweedie
  muz <- logit$linkinv(Xz %*% betaz + offz)
  mut <- tw$linkinv(Xt %*% betat + offt)
  mu <- as.numeric((1 - muz) * mut)  
  
  type <- match.arg(type)                                                            
  switch(type, response = mu, zero = muz, tweedie = mut)                                                            
})