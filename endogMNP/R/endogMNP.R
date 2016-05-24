endogMNP <- function (selForm, outForm, dataSet = parent.frame(), selBase = NULL, 
outBase = NULL, latent = FALSE, invcdf = FALSE, 
n.draws = 5000, p.var = "Inf", p.df = n.dim + 1, p.scale = 1, 
coef.start = 0, cov.start = 1, burnin = 0, thin = 0, verbose = FALSE, 
minConst = TRUE, trace = TRUE) 
{
    captureFormula <- function(formula, data = dataSet, latent = FALSE, 
							   invcdf = FALSE, n.draws = 5000, p.var = "Inf", 
							   p.df = n.dim + 1, p.scale = 1, coef.start = 0, cov.start = 1, 
							   burnin = 0, thin = 0, verbose = FALSE) {
        call <- match.call()
        mf <- match.call(expand = FALSE)
        mf$choiceX <- mf$cXnames <- mf$base <- mf$n.draws <- mf$latent <- mf$p.var <- 
		mf$p.df <- mf$p.scale <- mf$coef.start <- mf$invcdf <- mf$cov.start <- mf$verbose <- 
		mf$burnin <- mf$thin <- NULL
        mf[[1]] <- as.name("model.frame.default")
        mf$na.action <- "na.pass"
        mf <- eval.parent(mf)
        return(mf)
    }
    mf1 <- captureFormula(selForm, data = dataSet)
    mf2 <- captureFormula(outForm, data = dataSet)
    selY <- model.response(mf1)
    outY <- model.response(mf2)
    selX <- model.frame(mf1)[, -1]
    outX <- model.frame(mf2)[, -1]
    mf1$na.action <- "na.pass"
    Terms <- attr(mf1, "terms")
    selX <- model.matrix.default(Terms, mf1)
    mf2$na.action <- "na.pass"
    Terms <- attr(mf2, "terms")
    outX <- model.matrix.default(Terms, mf2)
    call <- match.call()
    tmp <- XandYmatrix.endogMNP(selX, selY, outX, outY, extra = TRUE, 
								verbose = verbose, base1 = selBase, base2 = outBase)
    Y <- tmp$Y
    lev1 <- tmp$lev1
    base1 <- tmp$base1
    p1 <- tmp$p1 - 1
    lev2 <- tmp$lev2
    base2 <- tmp$base2
    p2 <- tmp$p2 - 1
    selCovNum <- tmp$selCovNum
    outCovNum <- tmp$outCovNum
	notObsCat <- tmp$notObsCat
	whichNotObs <- tmp$whichNotObs
	if(notObsCat){
		n.dim <- p1 + p1*p2}
	else{
	n.dim <- p1 + (p1 + 1) * p2}
    if (verbose) 
	cat("\nThe base selection category is `", base1, "'.\n", 
		sep = "")
    if (verbose) 
	cat("\nThe base outcome category is `", base2, "'.\n", 
		sep = "")
    if (verbose) 
	cat("\nThe dimension of Sigma is ", n.dim, ".\n\n", sep = "")
    if (p1 < 1) 
	stop("The number of selection alternatives should be at least 2.")
    if (verbose) 
	cat("The total number of selection classes is ", p1 + 
		1, ".\n\n", sep = "")
    if (p2 < 1) 
	stop("The number of outcome alternatives should be at least 2.")
    if (verbose) 
	cat("The number of outcome classes is ", p2 + 1, 
		".\n\n", sep = "")
    X <- tmp$X
    
    selName <- tmp$selXnames
    outName <- tmp$outXnames
    coefnames <- NULL
    for (i in 2:(p1 + 1)) {
    	for (k in 1:length(selName)) {
            coefnames <- c(coefnames, paste("Sel-", selName[k], 
											" : ", lev1[i], sep = ""))
        }
    }
	if(notObsCat){
    for (l in 1:(p1 + 1)) {
		if(l != whichNotObs){
            for (i in 2:(p2 + 1)) {
            	for (k in 1:length(outName)) {
                coefnames <- c(coefnames, paste("Out-", outName[k], 
												" : ", lev2[i], "|", lev1[l], sep = ""))
            }}}}
	}
	else{
    for (l in 1:(p1 + 1)) {
    	for (i in 2:(p2 + 1)) {
			for (k in 1:length(outName)) {
					coefnames <- c(coefnames, paste("Out-", outName[k], 
												" : ", lev2[i], "|", lev1[l], sep = ""))
				}
				}
			}
				}
#number of covariates: 	
	if(notObsCat)
		n.cov <- p1 * selCovNum + (p1) * p2 * outCovNum
	else
	n.cov <- p1 * selCovNum + (p1 + 1) * p2 * outCovNum
## Cant handle NAs yet ## 	
    na.ind = 0
#number of observations 	
    n.obs <- nrow(X)/n.dim/n.cov
    if (verbose) {
        cat("The dimension of beta is ", n.cov, ".\n\n", sep = "")
        cat("The number of observations is ", n.obs, ".\n\n", 
            sep = "")
        if (sum(na.ind > 0) > 0) {
            if (sum(na.ind > 0) == 1) 
			cat("The observation ", (1:length(na.ind))[na.ind > 
				0], " is dropped due to missing values.\n\n", 
				sep = "")
            else {
                cat("The following ", sum(na.ind > 0), " observations are dropped due to missing values:\n", 
					sep = "")
                cat((1:length(na.ind))[na.ind > 0], "\n\n")
            }
        }
    }
    p.imp <- FALSE
    if (p.var == Inf) {
        p.imp <- TRUE
        p.prec <- diag(0, n.cov)
        if (verbose) 
		cat("Improper prior will be used for beta.\n\n")
    }
    else if (is.matrix(p.var)) {
        if (ncol(p.var) != n.cov || nrow(p.var) != n.cov) 
		stop("The dimension of `p.var' should be ", n.cov, 
			 " x ", n.cov, sep = "")
        if (sum(sign(eigen(p.var)$values) < 1) > 0) 
		stop("`p.var' must be positive definite.")
        p.prec <- solve(p.var)
    }
    else {
        p.var <- diag(p.var, n.cov)
        p.prec <- solve(p.var)
    }
    p.mean <- rep(0, n.cov)
    p.df <- eval(p.df)
    if (length(p.df) > 1) 
	stop("`p.df' must be a positive integer.")
    if (p.df < n.dim) 
	stop(paste("`p.df' must be at least ", n.dim, ".", sep = ""))
    if (abs(as.integer(p.df) - p.df) > 0) 
	stop("`p.df' must be a positive integer.")
    if (!minConst) {
        if (p.scale != 1) {
            p.scale <- 1
            warning("p.scale must be 1 when minConst=FALSE")
		}
		p.scale <- diag(rep(1,n.dim))
        }
    else {
        leadDimHold <- c(1, p1 + 1)
        for (i in 1:p1) {
            leadDimHold <- c(leadDimHold, p1 + 1 + (p2 * i))
        }
        if (!is.matrix(p.scale)) 
		p.scale <- diag(p.scale, n.dim)
        if (ncol(p.scale) != n.dim || nrow(p.scale) != n.dim) 
		stop("`p.scale' must be ", n.dim, " x ", n.dim, sep = "")
## form testerMat to see if p.scale is appropriately block-diagonal ## 		

		
        testerMat <- matrix(1, n.dim, n.dim)
        testerMat[1:p1, 1:p1] <- 0
		if(notObsCat){
		testerMat[(p1 + 1):n.dim, (p1 + 1):n.dim] <- 
		testerMat[(p1 + 1):n.dim, (p1 + 1):n.dim] - kronecker(diag(rep(1, (p1))), matrix(1, p2, p2))}
		else {
        testerMat[(p1 + 1):n.dim, (p1 + 1):n.dim] <- 
		testerMat[(p1 + 1):n.dim, (p1 + 1):n.dim] - kronecker(diag(rep(1, (p1 + 1))), matrix(1, p2, p2))}
        if (sum(abs(testerMat * p.scale)) != 0) {
            stop("`p.scale' must be block diagonal.  See manual.")
        }
 
        if (sum(sign(eigen(p.scale)$values) < 1) > 0) 
		stop("`p.scale' must be positive definite.")
        else if (trace == FALSE){
        	if (sum(abs(p.scale * testerMat)) != 0) {
            p.scale[leadDimHold, leadDimHold] <- 1
            warning("leading elements in the block diagonal p.scale will be set to 1.")
        } }
        else {
        	diagHld <- diag(p.scale)
        	if(sum(diagHld[1:p1]) != p1)
        	stop("The sub-matrices of `p.scale' must have trace equal to their dimension.")
        	diagHld <- diagHld[-(1:p1)]
        	while(length(diagHld) > 0){
        		if(sum(diagHld) != length(diagHld))
        		stop("The sub-matrices of `p.scale' must have trace equal to their dimension.")
        		diagHld <- diagHld[-(1:p2)]}
        	}
    }
    Signames <- NULL
    lev <- NULL
    for (i in 2:length(lev1)) {
        lev <- c(lev, paste("Sel-", lev1[i], sep = ""))
    }
	if(notObsCat){	
    for (i in 1:length(lev1)) {
        for (j in 2:length(lev2)) {
			if(i != whichNotObs){
            lev <- c(lev, paste("Out-", lev2[j], "|", lev1[i], 
								sep = ""))
			}
        }
    }}
	else{
		for (i in 1:length(lev1)) {
			for (j in 2:length(lev2)) {
				lev <- c(lev, paste("Out-", lev2[j], "|", lev1[i], 
									sep = ""))
			}
		}
	}
    for (j in 1:n.dim) {
        for (k in 1:n.dim) {
            if (j <= k) {
                Signames <- c(Signames, paste(lev[j], ":", lev[k], 
											  sep = ""))
            }
        }
    }
    if (length(coef.start) == 1) 
	coef.start <- rep(coef.start, n.cov)
    else if (length(coef.start) != n.cov) 
	stop(paste("The dimenstion of `coef.start' must be  ", 
			   n.cov, ".", sep = ""))
    if (!is.matrix(cov.start)) {
        cov.start <- diag(n.dim) * cov.start
        cov.start[1, 1] <- 1
    }
    else if (ncol(cov.start) != n.dim || nrow(cov.start) != n.dim) 
	stop("The dimension of `cov.start' must be ", n.dim, 
		 " x ", n.dim, sep = "")
    else if (sum(sign(eigen(cov.start)$values) < 1) > 0) 
	stop("`cov.start' must be a positive definite matrix.")
    else if (cov.start[1, 1] != 1) {
        cov.start[1, 1] <- 1
        warning("cov.start[1,1] will be set to 1.")
    }
    if (burnin < 0) 
	stop("`burnin' should be a non-negative integer.")
    if (thin < 0) 
	stop("`thin' should be a non-negative integer.")
    keep <- thin + 1
    if (latent) 
	n.par <- n.cov + n.dim * (n.dim + 1)/2 + n.dim * n.obs
    else n.par <- n.cov + n.dim * (n.dim + 1)/2
    if (verbose) 
	cat("Starting Gibbs sampler...\n")
## code NAs as -1 ##	
    Y[is.na(Y)] <- -1
	if(notObsCat){
    param <- .C("cMNPgibbsSel", as.integer(n.dim), as.integer(n.cov), 
				as.integer(p1), as.integer(p2), as.integer(n.obs), as.integer(n.draws), 
				as.double(p.mean), as.double(p.prec), as.integer(p.df), 
				as.integer(selCovNum), as.integer(outCovNum), as.double(p.scale),
				as.double(X), as.integer(Y), as.double(coef.start), 
				as.double(cov.start), as.integer(p.imp), as.integer(invcdf), 
				as.integer(burnin), as.integer(keep), as.integer(verbose), 
				as.integer(latent), as.integer(minConst), as.integer(trace),
				pdStore = double(n.par * floor((n.draws - burnin)/keep)), 
				PACKAGE = "endogMNP")$pdStore
	}
	else {
	
    param <- .C("cMNPgibbs", as.integer(n.dim), as.integer(n.cov), 
				as.integer(p1), as.integer(p2), as.integer(n.obs), as.integer(n.draws), 
				as.double(p.mean), as.double(p.prec), as.integer(p.df), 
				as.integer(selCovNum), as.integer(outCovNum), as.double(p.scale),
				as.double(X), as.integer(Y), as.double(coef.start), 
				as.double(cov.start), as.integer(p.imp), as.integer(invcdf), 
				as.integer(burnin), as.integer(keep), as.integer(verbose), 
				as.integer(latent), as.integer(minConst), as.integer(trace),
				pdStore = double(n.par * floor((n.draws - burnin)/keep)), 
				PACKAGE = "endogMNP")$pdStore
	}
    param <- matrix(param, ncol = n.par, nrow = floor((n.draws - 
													   burnin)/keep), byrow = TRUE)
    if (latent) {
        W <- array(as.vector(t(param[, (n.cov + n.dim * (n.dim + 1)/2 + 1):
							   (n.cov + n.dim * (n.dim + 1)/2 + n.dim * 
								n.obs)])), dim = c(n.dim, n.obs, floor((n.draws - 
								burnin)/keep)), dimnames = list(lev[-1], rownames(Y), NULL))
        param <- param[, 1:(n.cov + n.dim * (n.dim + 1)/2)]
    }
    else W <- NULL

    allNms <- c(coefnames, Signames)
    colnames(param) <- c(coefnames, Signames)
    Y[Y == -1] <- NA
    


        res <- list(call = call, param = param, x = X, y = Y, 
					n.dim = n.dim, n.obs = n.obs, 
					coefnames = coefnames, W = W, p.scale = p.scale, 
					n.cov = n.cov, nu0 = p.df, p.var = p.var, n.param = n.par, minConst = minConst, 
					n.dim1 = p1, n.dim2 = p2, n.rep = length(param[, 1]), 
					selForm = selForm, outForm = outForm, dataSet = dataSet, 
					selBase = base1, outBase = base2)
    
    class(res) <- "endogMNP"
    return(res)
}

