ordLORgee <-
function (formula = formula(data), data = parent.frame(), id = id, repeated = NULL, 
    link = "logit", bstart = NULL, LORstr = "category.exch", 
    LORem = "3way", LORterm = NULL, add = 0, homogeneous = TRUE, 
    restricted = FALSE, control = LORgee.control(), 
    ipfp.ctrl = ipfp.control(), IM = "solve") 
{
    options(contrasts=c("contr.treatment", "contr.poly"))
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
    LORstrs <- c("independence", "uniform", "category.exch", 
        "time.exch", "RC","fixed")
    icheck <- as.integer(match(LORstr, LORstrs, -1))
    if (icheck < 1) {
        stop("unknown local odds ratios structure")
    }
    if (LORstr == "independence" | LORstr == "fixed") {
        LORem <- NULL
    }
    else if (LORstr == "category.exch") {
        LORem <- "3way"
    }
    else {
        if (LORstr == "RC") 
            LORem <- "2way"
        if (LORem != "2way" & LORem != "3way") 
            stop("'LORem' must be '2way' or '3way'")
    }
    if (LORstr == "time.exch" | LORstr == "RC") {
        if (!is.logical(homogeneous)) 
            stop("'homogeneous' must be 'TRUE' or 'FALSE'")
        if (!is.logical(restricted)) 
            stop("'restricted' must be 'TRUE' or 'FALSE'")
        restricted <- if (!restricted) NULL else TRUE
    }
    else {
        homogeneous <- restricted <- NULL
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
    IMs <- c("cholesky", "solve","qr.solve")
    icheck <- as.integer(match(IM, IMs, -1))
    if (icheck < 1) 
        stop("unknown method for inverting a matrix")
    if (LORstr != "independence" & LORstr != "fixed") {
        data.model <- datacounts(Y, id, repeated, ncategories)
        marpars <- mmpar(LORem, LORstr, max(data.model$tp), homogeneous)
        LORem <- marpars$LORem
        LORstr <- marpars$LORstr
        LORterm <- fitmm(data.model, marpars,homogeneous, 
            restricted, add)
    }
    link <- as.character(link)
    links <- c("logit", "probit", "cloglog", "cauchit", "acl")
    icheck <- as.integer(match(link, links,-1))
    if (icheck < 1) { 
        stop("'link' must be \"logit\", \"probit\", \"cloglog\", \"cauchit\" or \"acl\"")
                    }
    if (is.null(bstart)) {
        if (link == "acl") family <- acat(reverse = TRUE, parallel = TRUE)
        if (link == "logit") family <- cumulative("logit", parallel = TRUE)
        if (link == "probit") family <- cumulative("probit", parallel = TRUE)
        if (link == "cloglog") family <- cumulative("cloglog", parallel = TRUE)
        if (link == "cauchit") family <- cumulative("cauchit", parallel = TRUE)
           mmcall <- match.call(expand.dots=FALSE)
           mmf <- match(c("formula", "data", "id", "repeated"), names(mmcall), 0L)
           mm <- mcall[c(1L, mmf)]
           mm$family <- family
           mm$control = vglm.control()
           mm[[1]] <- as.name("vglm")
           coeffs <- coef(eval(mm, parent.frame()))
           coeffs <- as.numeric(coeffs)
        if (!is.numeric(coeffs)) 
            stop("Please insert initial values")
        if (verbose) {
            cat("\nGEE FOR ORDINAL MULTINOMIAL RESPONSES\n")
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
    X_mat <- model.matrix(Terms, m)
    if(ncol(X_mat)>2){
    xnames <- colnames(X_mat)
    X_mat <- apply(X_mat[,-1],2,function(x) rep(x,each = ncategories-1))
    X_mat1 <- model.matrix(~factor(Intercept)-1)    
    X_mat <- cbind(X_mat1,X_mat)
    X_mat <- matrix(X_mat, ncol= ncol(X_mat), dimnames = NULL)
    xnames <- c(paste("beta0", 1:(ncategories - 1), sep = ""), xnames[-1])    
                      } else if(ncol(X_mat)==2) {
    xnames <- colnames(X_mat)
    X_mat <- rep(X_mat[,-1],each = ncategories-1)
    X_mat1 <- model.matrix(~factor(Intercept)-1)    
    X_mat <- cbind(X_mat1,X_mat)
    X_mat <- matrix(X_mat, ncol= ncol(X_mat), dimnames = NULL)
    xnames <- c(paste("beta0", 1:(ncategories - 1), sep = ""), xnames[-1]) 
                      } else {
    X_mat <- model.matrix(~factor(Intercept)-1)
    xnames <- c(paste("beta0", 1:(ncategories - 1), sep = ""))
                      }    
    if (link == "acl") {
        dummy <- ncategories - 1
        dummy.matrix <- diagmod(rep.int(1, dummy))
        dummy.matrix[upper.tri(dummy.matrix)] <- 1
        X_mat[, 1:dummy] <- kronecker(rep.int(1, nrow(X_mat)/dummy), 
            dummy.matrix)
        if (dummy != ncol(X_mat)) {
            X_mat[, -c(1:dummy)] <- X_mat[, -c(1:dummy)] * rep(dummy:1, 
                nrow(X_mat)/dummy)
        }
    }
    if (!is.null(bstart)) {
        coeffs <- as.numeric(bstart)
        if (length(coeffs) != ncol(X_mat)) 
            stop("Starting values and parameters vector differ in length")
        if (any(diff(coeffs[1:(ncategories - 1)]) < 0)) 
            stop("cutpoints are not increasing")
        if (verbose) {
            cat("\nGEE FOR ORDINAL MULTINOMIAL RESPONSES\n")
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
    fit$title <- "GEE FOR ORDINAL MULTINOMIAL RESPONSES"
    fit$version <- "version 1.5.1 modified 2015-03-09"
    fit$link <- if (link == "acl") 
        paste("Adjacent Category Logit")
    else paste("Cumulative", link, sep = " ")
    fit$local.odds.ratios <- list()
    fit$local.odds.ratios$structure <- LORstr
    fit$local.odds.ratios$model <- LORem
    fit$local.odds.ratios$homogeneous <- homogeneous
    fit$local.odds.ratios$restricted <- restricted
    fit$local.odds.ratios$theta <- fitmod$theta
    fit$terms <- Terms
    fit$contrasts <- attr(model.matrix(Terms, m), "contrasts")
    fit$convergence <- list()
    fit$convergence$niter <- fitmod$iter
    fit$convergence$criterion <- fitmod$crit[fitmod$iter]
    fit$convergence$conv <- fitmod$conv
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
    if (length(xnames) == (ncategories - 1)) 
        fit$pvalue <- NULL
    else {
        dummy <- 1:(ncategories - 1)
        waldts <- fit$coefficients[-dummy] %*% solve((fit$robust.variance)[-dummy, 
            -dummy])
        waldts <- waldts %*% fit$coefficients[-dummy]
        fit$pvalue <- 1 - pchisq(waldts, length(xnames) - length(dummy))
    }
    class(fit) <- "LORgee"
    fit
}