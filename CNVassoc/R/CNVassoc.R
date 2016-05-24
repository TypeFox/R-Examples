CNVassoc <- function (formula, data, subset, na.action, model = "multiplicative",
    family = "binomial", tol = 1e-06,
    max.iter = 30, emsteps = 0, verbose = FALSE, coef.start, sigma.start, alpha.start=1)
{
    cl <- match.call()
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "na.action"), names(mf), 0L) #offset
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    control <- unlist(lapply(mf, FUN = is.cnv))
    if (sum(control) == 0)
      stop("a variable of class 'cnv' should be included in the model")
    if (!is.null(attributes(mf)$na.action)){
      exclude <- as.integer(names(attributes(mf)$na.action))
      attr(mf[,control],"probabilities") <- attr(mf[,control],"probabilities")[-exclude,]
    }
    if (!missing(subset)){
      mforig <- mf
      attr.cnv <- attributes(mf[,control])
      mf <- match.call(expand.dots = FALSE)
      m <- match(c("formula", "data", "subset"), names(mf), 0L)
      mf <- mf[c(1L, m)]
      mf$drop.unused.levels <- TRUE
      mf[[1]] <- as.name("model.frame")
      mf <- eval(mf, parent.frame())
      attributes(mf[,control]) <- attr.cnv
      mf[,control] <- (mforig[,control])[rownames(mforig)%in%rownames(mf)]
    }
    mt <- attr(mf, "terms")
    model.type <- charmatch(model, c("multiplicative", "additive"))
    if (is.na(model.type))
        stop(" argument 'model' must be either 'multiplicative' or 'additive'")
    nIndiv <- NROW(mf)
    nameCNV <- names(mf)[control]
    mm <- model.matrix(mt, mf, contrasts)
    varPos <- which(colnames(mm) == nameCNV)
    var <- mf[, control]
    y <- mf[, 1]
    nVar <- NCOL(mf) - 2
    if (nVar != 0) {
        Xcov <- mm[, -c(1, varPos), drop = FALSE]
        nCov <- NCOL(Xcov)
    }
    else {
        Xcov <- NULL
        nCov <- 0
    }
    k <- attr(var, "k")
    w <- attr(var, "prob")
    num.copies <- attr(var, "num.copies")
    
    # arrange X & beta.start
    if (model.type == 1) {
        X <- array(dim = c(nIndiv, 1 + nCov, k), dimnames = list(1:nIndiv, c("CNVmult", colnames(Xcov)), 1:k))
        for (i in 1:k)
            X[, , i] <- cbind(rep(1, nIndiv), Xcov)
        if (nCov == 0)
            beta.start <- if (missing(coef.start)) rbind(rep(0, k)) else rbind(coef.start)
        else {
            beta.start <- matrix(c(rep(0, k), rep(c(0, rep(NA, k - 1)), nCov)), nrow = nCov + 1, ncol = k, byrow = TRUE)
            if (!missing(coef.start)) beta.start[1,]<-coef.start
        }        
    }
    else {
        X <- array(dim = c(nIndiv, 2 + nCov, k), dimnames = list(NULL, c("intercept", "CNVadd", colnames(Xcov)), 1:k))
        for (i in 1:k)
            X[, , i] <- cbind(rep(1, nIndiv), rep(num.copies[i], nIndiv), Xcov)
        beta.start <- matrix(0,nrow=nCov+2,ncol=k)
        beta.start[,-1]<-NA
        if (!missing(coef.start)) beta.start[,1]<-coef.start        
    }
    
    # perform CNVassoc 
    if (family=="binomial") {
      if (any(!y %in% c(0, 1)))
          stop("left hand side formula must be a variable with 0,1")
      if (emsteps > 0) {
          ans.temp <- EMlogistic(y, X, w, beta.start, tol = tol, max.iter = emsteps, verbose = verbose)
          beta.start <- ans.temp$beta
          variant <- ans.temp$variant
          ans <- NRlogistic(y, X, w, beta.start, variant, tol, max.iter = max.iter - emsteps, verbose)
      }
      else
          ans <- NRlogistic(y, X, w, beta.start, , tol, max.iter, verbose)
    } else {
      if(family=="gaussian") {
        if (missing(sigma.start)) 
            sigma.start <- sd(y)
        if (emsteps > 0) {
            ans.temp <- EMnorm(y, X, w, beta.start, sigma.start,
                tol = tol, max.iter = emsteps, verbose = verbose)
            beta.start <- ans.temp$beta
            sigma.start <- ans.temp$sigma
            variant <- ans.temp$variant
            ans <- NRnorm(y, X, w, beta.start, sigma.start, variant,
                tol = tol, max.iter = max.iter - emsteps, verbose = verbose)
        }
        else
            ans <- NRnorm(y, X, w, beta.start, sigma.start, , tol, max.iter, verbose)
      } else {
        if (family=="poisson") {
          if (any(y<0))
              stop("left hand side formula must be a non-negative integer")
          if (emsteps > 0) {
              ans.temp <- EMpoisson(y, X, w, beta.start, tol = tol, max.iter = emsteps, verbose = verbose)
              beta.start <- ans.temp$beta
              variant <- ans.temp$variant
              ans <- NRpoisson(y, X, w, beta.start, variant, tol, max.iter = max.iter - emsteps, verbose)
          }
          else
              ans <- NRpoisson(y, X, w, beta.start, , tol, max.iter, verbose)
        } else {
          if (family=="weibull") {
            if (!inherits(y, "Surv")) 
                stop("Response must be a survival object")
            cens <- y[ ,2]
            y <- y[ ,1]
            if (alpha.start<=0)
                stop("alpha.start must be a positive double")
            if (emsteps > 0) {
                ans.temp <- EMWeibull(y, cens, X, w, beta.start, alpha.start, tol = tol, max.iter = emsteps, verbose = verbose)
                beta.start <- ans.temp$beta
                alpha.start <- ans.temp$alpha
                variant <- ans.temp$variant
                ans <- NRWeibull(y, cens, X, w, beta.start, alpha.start, variant, tol, max.iter = max.iter - emsteps, verbose)
            }
            else
                ans <- NRWeibull(y, cens, X, w, beta.start, alpha.start, , tol, max.iter, verbose)          
          } else 
            stop("Family not implemented\n currently avalaible: 'gaussian', 'binomial', 'poisson' or 'weibull'")
        }
      }
    }

    varmat <- qr.solve(-ans$hessian)
    if (model.type == 1) {
        if (nCov == 0)
            vn <- paste("CNV", num.copies, sep = "")
        else vn <- c(paste("CNV", num.copies, sep = ""), rownames(ans$beta)[-1])
    }
    else {
        vn <- rownames(ans$beta)
    }
    if (family=="gaussian")
        vn <- c(vn, "sigma")
    if (family=="weibull")
        vn <- c(vn, "alpha")
    rownames(varmat) <- colnames(varmat) <- vn
    ans$call <- cl
    ans$y <- y
    if (family == "weibull")
      ans$cens <- cens
    ans$X <- X
    ans$w <- w
    ans$varmat <- varmat
    ans$formula <- formula
    ans$CNVname <- colnames(mf)[varPos]
    ans$CNV <- var
    if (missing(data))
        ans$data <- NULL
    else ans$data <- data
    ans$coefficients <- ans$beta
    ans$beta <- NULL
    colnames(ans$coefficients) <- paste("CNV", num.copies, sep = "")
    attr(ans, "family") <- family
    attr(ans, "model") <- model.type
    attr(ans, "nCov") <- nCov
    class(ans) <- "CNVassoc"
    ans

}

