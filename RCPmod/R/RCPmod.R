# This is package RCPmod 

"additive.logistic" <-
function(x) 
{
	tmp <- exp( x)
	tmp <- tmp / (1+sum( tmp))
	tmp <- c(tmp, 1-sum( tmp))
	
	return( tmp)

}


"AIC.regimix" <-
function (object, ..., k = 2) 
{
    habi.mod <- object
    p <- length(unlist(habi.mod$coefs))
    if (is.null(k)) 
        k <- 2
    star.ic <- -2 * habi.mod$logl + k * p
    return(star.ic)
}


"BIC.regimix" <-
function (object, ...) 
{
    habi.mod <- object
    p <- length(unlist(habi.mod$coefs))
    k <- log(habi.mod$n)
    star.ic <- -2 * habi.mod$logl + k * p
    return(star.ic)
}


"check.outcomes1" <-
function( outs) 
{
  nam <- colnames( outs)
  if( length( nam) == length( unique( nam)))
    return( length( nam))
  else
    return( FALSE)

}


"coef.regimix" <-
function (object, ...) 
{
    habi.mod <- object
    res <- list()
    res$alpha <- habi.mod$coefs$alpha
    res$tau <- matrix(habi.mod$coefs$tau, nrow = habi.mod$nRCP - 1, ncol = habi.mod$S)
    res$beta <- matrix(habi.mod$coefs$beta, nrow = habi.mod$nRCP - 1, ncol = habi.mod$p)
    return(res)
}


"extractAIC.regimix" <-
function (fit, scale = 1, k = 2, ...) 
{
    habi.mod <- fit
    n <- habi.mod$n
    edf <- length(unlist(coef(habi.mod)))
    if (is.null(k)) 
        k <- 2
    aic <- -2 * logLik(habi.mod) + k * edf
    return(c(edf, aic))
}


"globCIFinder" <-
function( x, en, alpha, nsim) 
{
	#this now works for both upper and lower CIs
	c <- uniroot( f=globErrorFn, interval=c(0.1,5), x=x, en=en, alpha=alpha, nsim=nsim)$root
		
	return( en*c)
	
}


"globErrorFn" <-
function( c1, x, en, alpha, nsim) 
{
	if( alpha > 0.5){
		tmp <- apply( x, 2, function(x) any( x-c1*en > 0))
		return( sum( tmp) / nsim - (1-alpha))
	}
	else{
		tmp <- apply( x, 2, function(x) any( x-c1*en < 0))
		return( sum( tmp) / nsim - alpha)
	}
}


"inv.logit" <-
function(x) 
{
	eta <- exp( x)
	mu <- eta / (1+eta)
	return(mu)
}


"logit" <-
function (x) 
	log( x / (1-x))


"logLik.regimix" <-
function (object, ...) 
{
    habi.mod <- object
    return(habi.mod$logl)
}


"margExpect" <-
function(pii, mu){
	wmu <- mu * pii
	mmu <- colSums( wmu)
	return( as.numeric( mmu))
}


"my.rmvnorm" <-
function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
    method = c("eigen", "svd", "chol")) 
{
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
        check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    sigma1 <- sigma
    dimnames(sigma1) <- NULL
    if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
        warning("sigma is numerically not symmetric")
    }
    method <- match.arg(method)
    if (method == "eigen") {
        ev <- eigen(sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
            warning("sigma is numerically not positive definite")
        }
        retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% 
            t(ev$vectors)
    }
    else if (method == "svd") {
        sigsvd <- svd(sigma)
        if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
            warning("sigma is numerically not positive definite")
        }
        retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    }
    else if (method == "chol") {
        retval <- chol(sigma, pivot = TRUE)
        o <- order(attr(retval, "pivot"))
        retval <- retval[, o]
    }
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
    retval <- sweep(retval, 2, mean, "+")
    colnames(retval) <- names(mean)
    retval
}


".onLoad" <-
function (libname, pkgname){
   # Generic DLL loader
   dll.path <- file.path( libname, pkgname, 'libs')
   if( nzchar( subarch <- .Platform$r_arch))
     dll.path <- file.path( dll.path, subarch)
   this.ext <- paste( sub( '.', '[.]', .Platform$dynlib.ext, fixed=TRUE), '$', sep='')

   dlls <- dir( dll.path, pattern=this.ext, full.names=FALSE)
   names( dlls) <- dlls
   if( length( dlls)) 
     lapply( dlls, function( x) library.dynam( sub( this.ext, '', x), package=pkgname, lib.loc=libname))
}


"orderFitted" <-
function( fm, simDat) 
{
#	require( gtools)
	
	RCPs <- attr( simDat, "RCP")
	posts <- fm$postProbs
	
	perms <- permutations( length( unique( RCPs)), length( unique( RCPs)))
	classErr <- rep( NA, ncol( perms))
	classErrRunnerUp <- classErr
	for( ii in 1:nrow( perms)){
		postsTMP <- posts[,perms[ii,]]
		postsTMP <- apply( postsTMP, 1, which.max)
		my.tab <- table( RCPs, postsTMP)
		classErr[ii] <- sum( diag( my.tab)) / sum( my.tab)
	}
	perms <- perms[which.max( classErr),]
	#coefs
	tau <- matrix( fm$coefs$tau, nrow=fm$nRCP-1, ncol=fm$S)
	tau <- rbind( tau, -colSums( tau))
	tau <- tau[perms,]
	beta <- matrix( fm$coefs$beta, nrow=fm$nRCP-1, ncol=fm$p)
	beta <- rbind( beta, 0)
	beta <- beta[perms,]
	beta <- beta - rep( beta[fm$nRCP,], each=3)
	fm$coefs <- list( alpha=fm$coefs$alpha, tau=as.numeric( tau[-fm$nRCP,]), beta=as.numeric( beta[-fm$nRCP,]))
	#scores
	fm$scores <- NULL
	#pis
	fm$pis <- fm$pis[,perms]
	#postProbs
	fm$postProbs <- fm$postProbs[,perms]
	#mus
	fm$mus <- fm$mus[perms,]
	#vcov
	fm$vcov <- NULL
	#order
	fm$perm <- perms
	#classification error
	fm$classErr <- max( classErr)
	fm$classErrRunnerUp <- max( classErr[-(which.max( classErr))])

	return( fm)

}


"orderPost" <-
function( new.fm=NULL, fm, RCPs=NULL, sample=NULL) 
{
#	require( gtools)
	G1 <- G2 <- NULL
	if( !is.null( new.fm))
		G <- G1 <- new.fm$nRCP
	if( !is.null( RCPs))
		G <- G2 <- length( unique( RCPs))
	if( sum( !is.null( c(G1,G2))) != 1){
		cat( "Problem with ordering -- provide new.fm *or* RCPs, but not both!\n")
		return( NULL)
	}
	perms <- permutations( G, G)
	
	if( !is.null( RCPs)){
		fm$postProbs <- matrix( 0, nrow=nrow( fm$postProbs), ncol=ncol( fm$postProbs))
		for( ii in 1:fm$nRCP)
			fm$postProbs[,ii] <- ifelse( RCPs==ii, 1, 0)
	}
	if( !is.null( sample))
		fm$postProbs <- fm$postProbs[sample,]
	classErr <- rep( NA, ncol( perms))
	for( ii in 1:nrow( perms)){
		my.tab <- t(fm$postProbs) %*% new.fm$postProbs[,perms[ii,]]
		classErr[ii] <- sum( diag( my.tab)) / sum( my.tab)
	}
	perms <- perms[which.max( classErr),]
	#coefs
	tau <- matrix( new.fm$coefs$tau, nrow=new.fm$nRCP-1, ncol=new.fm$S)
	tau <- rbind( tau, -colSums( tau))
	tau <- tau[perms,]
	beta <- matrix( new.fm$coefs$beta, nrow=new.fm$nRCP-1, ncol=new.fm$p)
	beta <- rbind( beta, 0)
	beta <- beta[perms,]
	beta <- beta - rep( beta[new.fm$nRCP,], each=3)
	new.fm$coefs <- list( alpha=new.fm$coefs$alpha, tau=as.numeric( tau[-new.fm$nRCP,]), beta=as.numeric( beta[-new.fm$nRCP,]))
	#scores
	new.fm$scores <- NULL
	#pis
	new.fm$pis <- new.fm$pis[,perms]
	#postProbs
	new.fm$postProbs <- new.fm$postProbs[,perms]
	#mus
	new.fm$mus <- new.fm$mus[perms,]
	#vcov
	new.fm$vcov <- NULL
	#order
	new.fm$perm <- perms
	#classification error
	new.fm$classErr <- max( classErr)
	new.fm$classErrRunnerUp <- max( classErr[-(which.max( classErr))])

	return( new.fm)



}


"plot.regimix" <-
function (x, ..., nsim = 100, alpha.conf = c(0.9, 0.95, 0.99)) 
{
    habi.mod <- x
    shad <- rev(seq(from = 0.8, to = 0.5, length = length(alpha.conf)))
    allResids <- matrix(NA, nrow = habi.mod$n, ncol = nsim)
    form <- habi.mod$titbit$form
    form[[2]] <- NULL
    form.X <- as.formula(form)
    X <- model.matrix(form.X, habi.mod$titbit$mf)
    alpha = habi.mod$coefs$alpha
    tau <- habi.mod$coefs$tau
    beta <- habi.mod$coefs$beta
    alpha.score <- as.numeric(rep(NA, habi.mod$S))
    tau.score <- as.numeric(matrix(NA, ncol = habi.mod$S, nrow = habi.mod$nRCP - 1))
    beta.score <- as.numeric(matrix(NA, ncol = ncol(X), nrow = habi.mod$nRCP - 1))
    scoreContri <- as.numeric(matrix(NA, ncol = length(unlist(habi.mod$coef)), nrow = habi.mod$n))
    pis <- as.numeric(matrix(NA, nrow = habi.mod$n, ncol = habi.mod$nRCP))
    mus <- as.numeric(matrix(NA, nrow = habi.mod$nRCP, ncol = habi.mod$S))
    logCondDens <- as.numeric(matrix(NA, nrow = habi.mod$n, ncol = habi.mod$nRCP))
    logls <- as.numeric(rep(NA, habi.mod$n))
    conv <- as.integer(0)
    pb <- txtProgressBar(min = 1, max = nsim, style = 3, char = "><(('> ")
    for (s in 1:nsim) {
        setTxtProgressBar(pb, s)
        newy <- as.matrix(simRCPdata(H = habi.mod$nRCP, 
            S = habi.mod$S, p = habi.mod$p, n = habi.mod$n, alpha = alpha, 
            tau = tau, beta = beta, X = X))
        tmp <- .Call("HABITAT_C", as.numeric(newy[, 1:habi.mod$S]), 
            as.numeric(X), as.integer(habi.mod$S), as.integer(habi.mod$nRCP), 
            as.integer(ncol(X)), as.integer(nrow(X)), alpha, 
            tau, beta, as.numeric(habi.mod$titbits$penalty[1]), 
            as.numeric(habi.mod$titbits$penalty[2]), alpha.score, 
            tau.score, beta.score, scoreContri, pis, mus, logCondDens, 
            logls, as.integer(20), as.integer(2), as.integer(1), 
            as.integer(1), as.numeric(0.01), as.numeric(0.01), 
            as.integer(conv), as.integer(4), as.numeric(TRUE), 
            PACKAGE = "RCPmod")
        allResids[, s] <- -2 * logls
    }
    cat("\n")
    allResidsSort <- apply(allResids, 2, sort)
    quants <- c(0.5, (1 - alpha.conf)/2, alpha.conf + (1 - alpha.conf)/2)
    envel <- t(apply(allResidsSort, 1, quantile, probs = quants, 
        na.rm = TRUE))
    sort.resid <- sort(habi.mod$residual)
    empQuant <- envel[, 1]
    diff <- sweep(envel[, -1], 1, empQuant, "-")
    realMeans <- (sort.resid + empQuant)/2
    realDiff <- sort.resid - empQuant
    par(mfrow = c(1, 2))
    plot(rep(realMeans, 1 + 2 * length(alpha.conf)), c(diff, 
        realDiff), sub = "Pointwise Confidence", 
        ylab = "Observed - Expected", xlab = "(Observed+Expected)/2", 
        type = "n")
    for (aa in length(alpha.conf):1) polygon(c(realMeans, rev(realMeans)), 
        c(diff[, aa], rev(diff[, aa + length(alpha.conf)])), 
        col = grey(shad[aa]), border = NA)
    points(realMeans, realDiff, pch = 20)
    abline(h = 0)
    globEnvel <- envel
    for (ii in 2:(length(alpha.conf) + 1)) globEnvel[, ii] <- globCIFinder(x = allResidsSort, en = envel[, ii], alpha = quants[ii], nsim = nsim)
    for (ii in 1 + (length(alpha.conf) + 1):(2 * length(alpha.conf))) globEnvel[, 
        ii] <- globCIFinder(x = allResidsSort, en = envel[, ii], 
        alpha = quants[ii], nsim = nsim)
    empQuant <- globEnvel[, 1]
    diff <- sweep(globEnvel[, -1], 1, empQuant, "-")
    realMeans <- (sort.resid + empQuant)/2
    realDiff <- sort.resid - empQuant
    plot(rep(realMeans, 1 + 2 * length(alpha.conf)), c(diff, 
        realDiff), sub = "Global Confidence", 
        ylab = "Observed - Expected", xlab = "(Observed+Expected)/2", 
        type = "n")
    for (aa in length(alpha.conf):1) polygon(c(realMeans, rev(realMeans)), 
        c(diff[, aa], rev(diff[, aa + length(alpha.conf)])), 
        col = grey(shad[aa]), border = NA)
    points(realMeans, realDiff, pch = 20)
    abline(h = 0)
    return(NULL)
}


"predict.regimix" <-
function (object, ..., newdata = NULL, nboot = 0, alpha = 0.95) 
{
    habi.mod <- object
    bootType <- "parametric"
    form <- habi.mod$titbit$form
    form.X <- as.formula(form)
    form.X[[2]] <- NULL
    mf <- habi.mod$titbits$mf
    if (is.null(newdata)) 
        mf.pred <- habi.mod$titbits$mf
    else mf.pred <- model.frame(form.X, data = as.data.frame(newdata))
    X <- model.matrix(form.X, mf.pred)
    S <- habi.mod$S
    G <- habi.mod$nRCP
    n <- nrow(X)
    p <- ncol(X)
    if (nboot > 0) 
        my.nboot <- nboot
    else my.nboot <- 0
    if (bootType == "parametric") 
        allCoBoot <- regibootParametric(fm = habi.mod, mf = mf, 
            nboot = my.nboot)
    if (bootType == "empirical") 
        allCoBoot <- regiboot(fm = habi.mod, nboot = my.nboot)
    if (is.null(allCoBoot)) 
        return(NULL)
    alphaBoot <- matrix(as.numeric(allCoBoot[, 1:S]), ncol = S)
    tauBoot <- as.numeric(allCoBoot[, S + 1:((G - 1) * S)], ncol = S)
    betaBoot <- as.numeric(allCoBoot[, S + (G - 1) * S + 1:((G - 
        1) * p)], ncol = S)
    alphaIn <- as.numeric(habi.mod$coefs$alpha)
    tauIn <- as.numeric(habi.mod$coef$tau)
    betaIn <- as.numeric(habi.mod$coef$beta)
    passType <- 1
    predCol <- G
    ptPreds <- as.numeric(matrix(NA, nrow = nrow(X), ncol = predCol))
    bootPreds <- as.numeric(array(NA, c(nrow(X), predCol, my.nboot)))
    conc <- as.numeric(NA)
    mysd <- as.numeric(NA)
    outcomes <- matrix(NA, nrow = nrow(X), ncol = S)
    tmp <- .Call("Habitat_predict_C", as.numeric(X), as.numeric(outcomes), 
        as.numeric(alphaIn), as.numeric(tauIn), as.numeric(betaIn), 
        as.numeric(alphaBoot), as.numeric(tauBoot), as.numeric(betaBoot), 
        as.integer(S), as.integer(G), as.integer(n), as.integer(p), 
        as.integer(nboot), ptPreds, bootPreds, as.integer(passType), 
        as.numeric(conc), as.numeric(mysd), PACKAGE = "RCPmod")
    ptPreds <- matrix(ptPreds, nrow = nrow(X), ncol = predCol)
    nam <- paste("RCP", 1:G, sep = "_")
    colnames(ptPreds) <- nam
    if (nboot > 0) {
        bootPreds <- matrix(bootPreds, nrow = nrow(X) * predCol, 
            ncol = my.nboot)
        bPreds <- list()
        row.exp <- rowMeans(bootPreds)
        tmp <- matrix(row.exp, nrow = nrow(X), ncol = predCol)
        bPreds$fit <- tmp
        tmp <- sweep(bootPreds, 1, row.exp, "-")
        tmp <- tmp^2
        tmp <- sqrt(rowSums(tmp)/(nboot - 1))
        tmp <- matrix(tmp, nrow = nrow(X), ncol = predCol)
        bPreds$ses <- tmp
        colnames(bPreds$fit) <- colnames(bPreds$ses) <- nam
        tmp <- apply(bootPreds, 1, quantile, probs = c(0, alpha) + 
            (1 - alpha)/2, na.rm = TRUE)
        tmp <- t(tmp)
        tmp <- array(tmp, c(nrow(X), predCol, 2), dimnames = list(NULL, 
            NULL, NULL))
        bPreds$cis <- tmp[, 1:predCol, ]
        dimnames(bPreds$cis) <- list(NULL, nam, c("lower", "upper"))
        ret <- list(ptPreds = ptPreds, bootPreds = bPreds$fit, 
            bootSEs = bPreds$ses, bootCIs = bPreds$cis)
    }
    else ret <- ptPreds
    gc()
    return(ret)
}


"print.regimix" <-
function (x, ...) 
{
    habi.mod <- x
    ret <- list()
    ret$Call <- habi.mod$call
    ret$coef <- coef(habi.mod)
    return(ret)
}


"regiboot" <-
function (fm, nboot) 
{
    if (nboot > 0) {
        boot.estis <- matrix(NA, nrow = nboot, ncol = length(unlist(fm$coef)))
        fm$titbits$control$boot <- TRUE
        fm$titbits$control$reltol <- max(1e-05, fm$titbits$control$reltol)
        pb <- txtProgressBar(min = 1, max = nboot, style = 3, 
            char = "><(('> ")
        for (ii in 1:nboot) {
            setTxtProgressBar(pb, ii)
            samp <- sample(1:fm$n, fm$n, replace = TRUE)
            samp.data <- fm$titbits$mf[samp, ]
            dumbOut <- capture.output(samp.fm <- regimix(form = fm$titbits$form, 
                data = samp.data, nRCP = fm$nRCP, control = fm$titbits$control, 
                inits = unlist(fm$coef)))
            boot.estis[ii, ] <- unlist(samp.fm$coefs)
        }
        cat("\n")
        fm$titbits$control$boot <- NULL
        return(boot.estis)
    }
    else {
        boot.estis <- matrix(unlist(fm$coef), nrow = 1)
    }
}


"regibootParametric" <-
function( fm, mf, nboot) 
{
	if( nboot > 0){
		S <- fm$S
		G <- fm$nRCP
		n <- fm$n
		p <- fm$p
		pstar <- length( unlist( fm$coef))
		
		betaID <- S+(G-1)*S + 1:((G-1)*p)
		if( is.null( fm$vcov)){
			cat( "An estimate of the variance matrix for regression parameters is required. Please run fm$vcov <- vcov(), see ?vcov.regimix for help\n")
			return( NULL)
		}
		
		allCoBoot <- matrix( rep( unlist( fm$coefs), each=nboot), nrow=nboot, ncol=pstar)
		allCoBoot[,betaID] <- my.rmvnorm( n=nboot, mean=as.numeric( fm$coefs$beta), sigma=fm$vcov[betaID,betaID], method="eigen")
		
		return( allCoBoot)
	}
	else{
		boot.estis <- matrix( unlist( fm$coef), nrow=1)
		return( boot.estis)
	}
}


"regimix" <-
function (form = NULL, data, nRCP = 3, control = list(), inits = "random") 
{
    call <- match.call()
    form <- as.formula(form)
    data <- as.data.frame(data)
    G <- as.integer(nRCP)
    if (is.null(control$boot)) 
        mf <- model.frame(form, data = data, na.action = na.omit)
    else mf <- data
    outcomes <- model.response(mf)
    form.X <- form
    form.X[[2]] <- NULL
    form.X <- as.formula(form.X)
    X <- model.matrix(form.X, mf)
    S <- check.outcomes1(outcomes)
    if (!S) {
        print("Two species have the same name -- exitting now")
        return(NULL)
    }
    cat("There are", S, "species", "\n")
    n.tot <- nrow(data)
    na.vec <- !apply(cbind(outcomes, X), 1, function(x) any(is.na(x)))
    n <- sum(na.vec)
    cat("There are", n, "fully present observations and", n.tot, 
        "observations in total", "\n")
    outcomes <- outcomes[na.vec, ]
    X <- X[na.vec, ]
    if (!("maxit" %in% names(control))) 
        control$maxit <- 500
    if (!("maxitInner" %in% names(control))) 
        control$maxitInner <- 100
    if (!("trace" %in% names(control))) 
        control$trace <- 1
    if (!("nreport" %in% names(control))) 
        control$nreport <- 1
    if (!("abstol" %in% names(control))) 
        control$abstol <- 1e-05
    if (!("reltol" %in% names(control))) 
        control$reltol <- sqrt(.Machine$double.eps)
    if (!("loglOnly" %in% names(control))) 
        loglOnly <- FALSE
    if (!("penalty" %in% names(control))) 
        penalty <- 0
    else {
        if (control$penalty < 0) {
            cat("Supplied penalty is negative, reverting to the default\n")
            penalty <- 0
        }
        else penalty <- control$penalty
    }
    if (!("maxGSiters" %in% names(control))) 
        GSiters <- 0
    else GSiters <- control$maxGSiters

    alpha <- beta <- tau <- -99999
    if (inits[1] == "hclust" | inits[1]=="random") {
        tmp <- dist(outcomes, method = "manhattan")
        tmp1 <- hclust(tmp, method = "ward")
        tmpGrp <- cutree(tmp1, G)
        mus <- apply(outcomes, 2, function(x) tapply(x, tmpGrp, mean))
        mus <- ifelse(mus <= 0.1, 0.1, mus)
        mus <- ifelse(mus >= 0.9, 0.9, mus)
        logitMus <- logit(mus)
        alpha <- as.numeric(colMeans(logitMus))
        tau <- as.numeric((logitMus - rep(alpha, each = G))[1:(G - 1), ])
        beta <- as.numeric(matrix(0, ncol = ncol(X), nrow = G - 1))
    }
    if (inits[1] == "hclust")
        cat( "Obtaining initial values for species' model from simple clustering algorithm", "\n")
    if (inits[1] == "random") {
        my.sd <- 0.1
        alpha <- alpha + rnorm(S, sd = my.sd)
        tau <- tau + as.numeric(matrix(rnorm((G - 1) * S, sd = my.sd), ncol = G - 1))
        beta <- beta + as.numeric(matrix(rnorm((G - 1) * ncol(X), mean = 0, sd = my.sd), ncol = ncol(X), nrow = G - 1))
    }
    if( alpha == -99999) {
        cat("Using supplied initial values (unchecked). Responsibility is entirely the users!\n")
        alpha <- inits[1:S]
        tau <- inits[S + 1:((G - 1) * S)]
        beta <- inits[S + (G - 1) * S + 1:((G - 1) * ncol(X))]
    }
    inits <- c(alpha, tau, beta)
    alpha.score <- as.numeric(rep(NA, S))
    tau.score <- as.numeric(matrix(NA, ncol = S, nrow = G - 1))
    beta.score <- as.numeric(matrix(NA, ncol = ncol(X), nrow = G - 
        1))
    scoreContri <- as.numeric(matrix(NA, ncol = length(inits), 
        nrow = n))
    pis <- as.numeric(matrix(NA, nrow = n, ncol = G))
    mus <- as.numeric(matrix(NA, nrow = G, ncol = S))
    logCondDens <- as.numeric(matrix(NA, nrow = n, ncol = G))
    logls <- as.numeric(rep(NA, n))
    conv <- as.integer(0)
    tmp <- .Call("HABITAT_C", as.numeric(outcomes), as.numeric(X), 
        as.integer(S), as.integer(G), as.integer(ncol(X)), as.integer(nrow(X)), 
        alpha, tau, beta, as.numeric(penalty), as.numeric(-99999), 
        alpha.score, tau.score, beta.score, scoreContri, pis, 
        mus, logCondDens, logls, as.integer(control$maxit), as.integer(control$maxitInner), 
        as.integer(control$trace), as.integer(control$nreport), 
        as.numeric(control$abstol), as.numeric(control$reltol), 
        as.integer(conv), as.integer(GSiters), as.numeric(loglOnly), PACKAGE = "RCPmod")
    coefs <- list(alpha = alpha, tau = tau, beta = beta)
    scores <- list(alpha = alpha.score, tau = tau.score, beta = beta.score)
    scoreContri <- matrix(scoreContri, nrow = n, ncol = length(inits))
    pis <- matrix(pis, ncol = G)
    mus <- matrix(mus, nrow = G)
    logCondDens <- matrix(logCondDens, ncol = G)
    postProbs <- pis * exp(logCondDens)
    postProbs <- postProbs/rowSums(postProbs)
    residuals <- -2 * logls
    k <- length(unlist(coefs))
    bic <- -2 * tmp + log(n) * k
    aic <- -2 * tmp + 2 * k
    entro <- postProbs * log(postProbs)
    EN <- -sum(entro)
    ICL <- bic + 2 * EN
    ret <- list(logl = tmp, coefs = coefs, scores = scores, pis = pis, 
        postProbs = postProbs, mus = mus, residuals = residuals, 
        n = n, S = S, p = ncol(X), nRCP = nRCP, BIC = bic, 
        AIC = aic, ICL = ICL, start.vals = inits, call = call, 
        titbits = list(mf = mf, form = form, penalty = penalty, 
            control = control, scoreContri = scoreContri))
    gc()
    class(ret) <- "regimix"
    return(ret)
}


"residuals.regimix" <-
function( x) 
{
	resid <- x$residuals
	return( resid)
}


"simRCPdata" <-
function (H = 3, S = 20, p = 3, n = 200, alpha = NULL, tau = NULL, 
    beta = NULL, X = NULL) 
{
    if (is.null(alpha) | length(alpha) != S) {
        cat("Scenario #1 values for alpha (identically 0)\n")
        alpha <- rep(0, S)
    }
    if (is.null(tau) | length(tau) != (H - 1) * S) {
        cat("Scenario #1 values for tau (+/- log(3))\n")
        tau <- sample(x = c(-1, 1), size = (H - 1) * S, replace = TRUE)
        tau <- tau * log(3)
    }
    tau <- matrix(as.numeric(tau), nrow = H - 1)
    if (is.null(beta) | length(beta) != (H - 1) * p) {
        cat("Scenario #1 values for beta\n")
        beta <- as.numeric(c(0, 0, 0.4, 0, -0.2, 1))
        beta <- beta
    }
    beta <- matrix(as.numeric(beta), nrow = H - 1)
    sppNames <- paste("spp", 1:S, sep = "")
    if (is.null(X)) {
        cat("creating a design matrix with random numbers\n")
        X <- cbind(1, matrix(runif(n * (p - 1), min = -10, max = 10), 
            nrow = n))
        colnames(X) <- c("intercept", paste("x", 1:(p - 1), sep = ""))
    }
    etaPi <- X %*% t(beta)
    pis <- t(apply(etaPi, 1, additive.logistic))
    habis <- apply(pis, 1, function(x) sample(1:H, 1, FALSE, 
        x))
    tau <- rbind(tau, -colSums(tau))
    etaMu <- tau + rep(alpha, each = H)
    mu <- inv.logit(etaMu)
    fitted <- mu[habis, ]
    outcomes <- matrix(rbinom(n * S, 1, fitted), nrow = n, ncol = S)
    colnames(outcomes) <- paste("spp", 1:S, sep = "")
    res <- as.data.frame(cbind(outcomes, X))
    attr(res, "RCPs") <- habis
    attr(res, "pis") <- pis
    attr(res, "alpha") <- alpha
    attr(res, "tau") <- tau[-H, ]
    attr(res, "beta") <- beta
    attr(res, "mu") <- mu
    return(res)
}


"summary.regimix" <-
function (object, ...) 
{
    habi.mod <- object
    if (is.null(habi.mod$vcov)) {
        habi.mod$vcov <- matrix(NA, nrow = length(unlist(habi.mod$coef)), 
            ncol = length(unlist(habi.mod$coef)))
        cat("No variance matrix has been supplied\n")
    }
    cat("Standard errors for alpha and tau parameters may be (are likely to be) misleading\n")
    res <- cbind(unlist(habi.mod$coefs), sqrt(diag(habi.mod$vcov)))
    res <- cbind(res, res[, 1]/res[, 2])
    res <- cbind(res, 2 * (1 - pnorm(abs(res[, 3]))))
    colnames(res) <- c("Estimate", "SE", "z-score", "p")
    return(res)
}


"vcov.regimix" <-
function (object, ..., method = "simple", nboot = 100) 
{
    habi.mod <- object
    if (!method %in% c("EmpInfo", "simple", "Richardson", "boot")) {
        cat("Unknown method to calculate variance matrix, options are: 'EmpInfo', 'simple', 'Richardson' and 'boot'\n")
        return(NULL)
    }
    if (method == "boot" & nboot <= 2) {
        cat("Bootstrap with (less than) 2 resamples is unlikely to work -- even if it does it will produce highly variable results.\n")
        return(NULL)
    }
    fm <- habi.mod
    outcomes <- model.response(fm$titbits$mf)
    form.X <- fm$titbits$form
    form.X[[2]] <- NULL
    form.X <- as.formula(form.X)
    X <- model.matrix(form.X, fm$titbits$mf)
    if (method %in% c("simple", "Richardson")) {
 #       require(numDeriv)
        pis <- matrix(as.numeric(NA), nrow = nrow(fm$pis), ncol = ncol(fm$pis))
        mus <- matrix(as.numeric(NA), nrow = nrow(fm$mus), ncol = ncol(fm$mus))
        logCondDens <- matrix(as.numeric(NA), nrow = nrow(fm$postProbs), 
            ncol = ncol(fm$postProbs))
        logls <- as.numeric(rep(NA, length(fm$residuals)))
        my.fun <- function(x) {
            alpha <- x[1:fm$S]
            tau <- x[fm$S + 1:((fm$nRCP - 1) * fm$S)]
            beta <- x[fm$S + (fm$nRCP - 1) * fm$S + 1:((fm$nRCP - 1) * ncol(X))]
            alpha.score <- as.numeric(rep(NA, fm$S))
            tau.score <- as.numeric(matrix(NA, nrow = fm$nRCP - 1, ncol = fm$S))
            beta.score <- as.numeric(matrix(NA, nrow = fm$nRCP - 1, ncol = ncol(X)))
            hess <- as.numeric(matrix(NA, length(x), length(x)))
            tmp <- .Call("HABITAT_C", as.numeric(outcomes), as.numeric(X), 
                as.integer(fm$S), as.integer(fm$nRCP), as.integer(ncol(X)), 
                as.integer(nrow(X)), alpha, tau, beta, as.numeric(fm$titbits$penalty[1]), 
                as.numeric(100), alpha.score, tau.score, beta.score, 
                hess, pis, mus, logCondDens, logls, as.integer(1000), 
                as.integer(0), as.integer(0), as.integer(0), 
                as.numeric(0.01), as.numeric(0.01), as.integer(1), 
                as.integer(0), as.numeric(TRUE), PACKAGE = "RCPmod")
            tmp1 <- c(alpha.score, tau.score, beta.score)
            return(tmp1)
        }
#        require(numDeriv)
        if (method == "Richardson") 
            cat("Method is Richardson approximation, may take a while if there is lots of data, species and/or RCP types\n")
        hess <- jacobian(my.fun, unlist(fm$coefs), method = method)
        hess <- 0.5 * (hess + t(hess))
        vcov.mat <- -solve(hess)
    }
    if (method == "boot") {
        boot.estis <- regiboot(fm = fm, nboot = nboot)
        vcov.mat <- var(boot.estis)
    }
    if (method == "EmpInfo") {
        p <- length(unlist(fm$coef))
        vcov.mat <- matrix(NA, nrow = p, ncol = p)
        betaID <- fm$S + (fm$nRCP - 1) * fm$S + 1:((fm$nRCP - 1) * fm$p)
        betaScoreContri <- fm$titbits$scoreContri[, fm$S + (fm$nRCP - 1) * fm$S + 1:((fm$nRCP - 1) * ncol(X))]
        beta.score <- fm$scores$beta
        tmp <- t(apply(betaScoreContri, 1, function(x) x %o% 
            x))
        tmp <- colSums(tmp)
        tmp <- matrix(tmp, nrow = length(beta.score), ncol = length(beta.score))
        tmp <- tmp - unlist(beta.score) %o% unlist(beta.score)/fm$n
        betaVcov <- solve(tmp)
        vcov.mat[betaID, betaID] <- betaVcov
    }
    rm(fm, habi.mod)
    return(vcov.mat)
}

