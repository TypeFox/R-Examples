nomLORgee <-
function (formula = formula(data), data = parent.frame(), id = id, repeated = NULL, 
    bstart = NULL, LORstr = "time.exch", LORem = "3way",
    LORterm = NULL, add = 0, homogeneous = TRUE, 
    control = LORgee.control(), ipfp.ctrl = ipfp.control(), IM = "solve") 
{
   options(contrasts=c("contr.treatment", "contr.poly"))
   restricted <- NULL
   link <- "bcl"
   call <- match.call() 
   mcall <- match.call(expand.dots=FALSE)
   mf <- match(c("formula", "data", "id", "repeated"), names(mcall), 0L)
   m <- mcall[c(1L, mf)]
   if (is.null(m$id)) 
   m$id <- as.name("id")
   m[[1]] <- as.name("model.frame")
   m <- eval(m, envir = parent.frame())
   Terms <- attr(m, "terms") 
    if(attr(Terms,"intercept")!=1) 
       stop("an intercept must be included")
   Y <- as.numeric(factor(model.response(m)))
    if (is.null(Y)) {
        stop("response variable not found")
    }
    ncategories <- nlevels(factor(Y))
    if (ncategories <= 2) 
        stop("The response variable should have more than 2 categories")
    id <- model.extract(m, "id")
    if (is.null(id)) {
        stop("'id' variable not found")
    }
    if (length(id) != length(Y)) 
        stop("response variable and 'id' are not of same length")
    repeated <- model.extract(m, "repeated")
    if (is.null(repeated)) {
      index <- order(unlist(split(1:length(id),id)))
      repeated <- c(unlist(sapply(unlist(lapply(split(id, id), length)), function(x) 1:x)))
      repeated <- repeated[index]
    }
    if (length(repeated) != length(Y)) 
        stop("response variable and 'repeated' are not of same length")
    id <- as.numeric(factor(id))
    repeated <- as.numeric(factor(repeated))
    if(all(id==repeated)) 
         stop("'repeated' and 'id' must not be equal")
    dummy <- split(repeated, id)
    if (any(unlist(lapply(dummy, length)) != unlist(lapply(lapply(dummy, 
        unique), length)))) 
        stop("'repeated' does not have unique values per 'id'")
    offset <- model.extract(m, "offset")
    if (length(offset) <= 1) 
    offset <- rep(0, length(Y))
    if (length(offset) != length(Y)) 
        stop("response variable and 'offset' are not of same length")
    offset <- as.double(offset)
    LORstrs <- c("independence", "time.exch", "RC", "fixed")
    icheck <- as.integer(match(LORstr, LORstrs, -1))
    if (icheck < 1) {
        stop("unknown local odds ratio structure")
    }
    if (LORstr == "independence" | LORstr == "fixed") {
        LORem <- NULL
    } else {
        if (LORstr == "RC") 
            LORem <- "2way"
        if (LORem != "2way" & LORem != "3way") 
            stop("'LORem' must be '2way' or '3way'")
    }
    if (LORstr == "time.exch" | LORstr == "RC") {
        if (!is.logical(homogeneous)) 
            stop("'homogeneous' must be 'TRUE' or 'FALSE'")
    }
    else {
        homogeneous <- NULL
    }
    if (LORstr == "independence" | LORstr == "fixed") {
        add <- NULL
    }
    else {
        if (!is.numeric(add) | add < 0) 
            stop("'add' must be >=0")
    }
    ipfp.ctrl <- ipfp.ctrl
    control <- control
    verbose <- control$verbose
    IMs <- c("cholesky", "solve", "qr.solve" )
    icheck <- as.integer(match(IM, IMs, -1))
    if (icheck < 1) 
        stop("unknown method for inverting a matrix")
    if (LORstr != "independence" & LORstr != "fixed") {
        data.model <- datacounts(Y, id, repeated, ncategories)
        marpars <- mmpar(LORem, LORstr, max(data.model$tp), homogeneous)
        LORem <- marpars$LORem
        LORstr <- marpars$LORstr
        LORterm <- fitmm(data.model, marpars, homogeneous, 
            NULL, add)
    }    
    if (is.null(bstart)) {
        family <- multinomial(refLevel = ncategories)
        mmcall <- match.call(expand.dots=FALSE)
        mmf <- match(c("formula", "data", "id", "repeated"), names(mmcall), 0L)
        mm <- mcall[c(1L, mmf)]
        mm$family <- family
        mm$control = vglm.control()
        mm[[1]] <- as.name("vglm")
        coeffs <- coef(eval(mm, parent.frame()))
        coeffs <- c(matrix(coeffs, ncol = ncategories - 
            1, byrow = TRUE))
        coeffs <- as.numeric(coeffs)
        if (!is.numeric(coeffs)) 
            stop("Please insert initial values")
        if (verbose) {
            cat("\nGEE FOR NOMINAL MULTINOMIAL RESPONSES\n")
            cat("\nrunning 'vglm' function to get initial regression estimates\n")
            print(matrix(coeffs, ncol = 1, dimnames = list(1:length(coeffs), 
               "Initial.Values")))
        }
    }
    Y <- rep(Y, each=ncategories-1)
    Intercept <- rep.int(seq(ncategories - 1),length(id))
    Y <- as.numeric(Y == Intercept)
    id <- rep(id, each = ncategories - 1)
    repeated <- rep(repeated, each = ncategories - 1)    
    offset <- rep(offset, each = ncategories-1)
    Xinit_mat <- model.matrix(Terms, m)
    xxnames <- colnames(Xinit_mat)
    Xinit_mat <- apply(Xinit_mat,2,function(x) rep(x,each=ncategories-1))
    X_mat <- model.matrix(~factor(Intercept)-1)
    if(ncol(Xinit_mat)>1) {
      Xinit_mat <- cbind(X_mat,Xinit_mat[,-1]) } else {
      Xinit_mat <- X_mat
                          }   
    if (ncol(Xinit_mat) != (ncategories - 1)) {
        X_inter <- X_mat
        for (i in ncategories:ncol(Xinit_mat)) X_mat <- cbind(X_mat, 
            X_inter * Xinit_mat[, i])
    }
    X_mat <- matrix(X_mat, ncol = ncol(X_mat), dimnames = NULL)
    X_mat <- X_mat[, c(matrix(seq(ncol(X_mat)), ncol = ncategories - 
        1, byrow = TRUE))]    
    if (!is.null(bstart)) {
        coeffs <- as.numeric(bstart)
        if (length(coeffs) != ncol(X_mat)) 
            stop("'bstart' and 'beta' differ in length")
        if (verbose) {
            cat("\nGEE FOR NOMINAL MULTINOMIAL RESPONSES\n")
            cat("\nuser's initial regression estimate\n")
            print(matrix(coeffs, ncol = 1, dimnames = list(1:length(coeffs), 
                "Initial.Values")))
        }
    }
    ordindex <- order(id,repeated)
    Y <- Y[ordindex]
    X_mat <- X_mat[ordindex,]
    id <- id[ordindex]
    repeated <- repeated[ordindex]
    offset <- offset[ordindex]
    fitmod <- fitLORgee(Y, X_mat, coeffs, ncategories, id, repeated, offset, 
    link, LORterm, marpars, ipfp.ctrl, control, IM, LORem = LORem, 
    LORstr = LORstr, add)
    fit <- list()
    fit$call <- call
    fit$title <- "GEE FOR NOMINAL MULTINOMIAL RESPONSES"
    fit$version <- "version 1.5.1 modified 2015-03-09"
    fit$link <- c("Baseline Category Logit")
    fit$local.odds.ratios <- list()
    fit$local.odds.ratios$structure <- LORstr
    fit$local.odds.ratios$model <- LORem
    fit$local.odds.ratios$homogeneous <- homogeneous
    fit$local.odds.ratios$theta <- fitmod$theta
    fit$terms <- Terms
    fit$contrasts <- attr(model.matrix(Terms, m), "contrasts")
    fit$convergence <- list()
    fit$convergence$niter <- fitmod$iter
    fit$convergence$criterion <- fitmod$crit[fitmod$iter]
    fit$convergence$conv <- fitmod$conv
    xnames <- paste("beta0", 1:(ncategories - 
           1), sep = "") 
    if(length(xxnames)>1) {
        xxnames <- c(xnames,xxnames[-1]) 
     if (length(xxnames) > length(xnames)) {
        for (i in 1:((ncol(X_mat) - ncategories + 1)/(ncategories - 
            1))) xnames <- c(xnames, paste(xxnames[i + ncategories - 
            1], 1:(ncategories - 1), sep = ":"))
                                           }
    xnames <- xnames[c(matrix(seq(ncol(X_mat)), ncol = ncategories - 
        1, byrow = TRUE))]
                          } 
    fit$coefficients <- fitmod$beta_mat[, fitmod$iter + 1]
    names(fit$coefficients) <- xnames
    fit$linear.predictors <- matrix(fitmod$linear.predictor, 
        ncol = ncategories - 1, byrow = TRUE)
    rownames(fit$linear.predictors) <- 1:nrow(fit$linear.predictors)
    colnames(fit$linear.predictors) <- 1:(ncategories - 1)
    fitted.values <- fitmod$fitted.values
    fitted.values.1 <- matrix(fitted.values, ncol = ncategories - 
        1, byrow = TRUE)
    fitted.values.2 <- 1 - rowSums(fitted.values.1)
    fitted.values <- cbind(fitted.values.1, fitted.values.2)
    rownames(fitted.values) <- 1:nrow(fitted.values.1)
    colnames(fitted.values) <- 1:ncategories
    fit$fitted.values <- fitted.values
    fit$residuals <- matrix(fitmod$residuals, ncol = ncategories - 
        1, byrow = TRUE)
    rownames(fit$residuals) <- 1:nrow(fit$residuals)
    colnames(fit$residuals) <- 1:(ncategories - 1)
    y <- Y
    y <- apply(matrix(y, ncol = ncategories - 1, byrow = TRUE), 
        1, function(x) which(x == 1))
    y <- as.numeric(y)
    y[is.na(y)] <- ncategories
    fit$y <- y
    fit$nobs <- length(y)
    fit$max.id <- max(unique(id))
    fit$clusz <- unlist(lapply(split(id, id), length))/(ncategories - 
        1)
    fit$id <- rep(1:fit$max.id,as.numeric(fit$clusz))
    fit$robust.variance <- fitmod$robust
    dimnames(fit$robust.variance) <- list(xnames, xnames)
    fit$naive.variance <- fitmod$naive
    dimnames(fit$naive.variance) <- list(xnames, xnames)
    fit$xnames <- xnames
    fit$categories <- ncategories
    fit$occasions <- sort(unique(repeated))
    fit$LORgee.control <- control
    fit$ipfp.control <- ipfp.ctrl
    fit$inverse.method <- IM
    fit$adding.constant <- add
    if (control$TRACE) {
        fit$trace <- list()
        fit$trace$coeffs <- fitmod$beta_mat
        fit$trace$crit <- fitmod$crit
    }
    if (length(xxnames) == (ncategories - 1)) 
        fit$pvalue <- NULL
    else {
        dummy <- seq(1, length(xxnames), ncategories - 1)
        waldts <- fit$coefficients[-dummy] %*% solve((fit$robust.variance)[-dummy, 
            -dummy])
        waldts <- waldts %*% fit$coefficients[-dummy]
        fit$pvalue <- 1 - pchisq(waldts, length(xxnames) - length(dummy))
    }
    class(fit) <- "LORgee"
    fit
}