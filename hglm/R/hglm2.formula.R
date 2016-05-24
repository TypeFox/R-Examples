`hglm2.formula` <-
	function(meanmodel, data = NULL, family = gaussian(link = identity),
             rand.family = gaussian(link = identity), method = "EQL",
             conv = 1e-6, maxit = 50, startval = NULL, X.disp = NULL, disp = NULL,
             link.disp = "log", weights = NULL, fix.disp = NULL, offset = NULL, 
             sparse = TRUE, vcovmat = FALSE, calc.like = FALSE, RandC = NULL,
			 bigRR = FALSE, verbose = FALSE, ...) {
         
Call <- match.call(expand.dots = FALSE)

### check the formulae ###
if (!inherits(meanmodel, "formula") || length(meanmodel) < 3) stop("\n Mean model must be a formula of the form \"response ~ fixd + (random)\"")
MainResponse <- all.vars(meanmodel)[1]
MainTerms <- terms(meanmodel)
Intercept <- attributes(MainTerms)$intercept
MainTerms <- attributes(MainTerms)$term.labels
RandTerms <- grepl("\\|",MainTerms)
NrRef <- sum(RandTerms)
if (NrRef == 0) stop("meanmodel must contain one or more random effects")
if (length(MainTerms) == NrRef & Intercept == 0) stop("Model must contain at least one fixed effect, e.g. an intercept")
MainTerms <- c(Intercept, MainTerms)
RandTerms <- c(FALSE, RandTerms)
fixed <- as.formula(paste(MainResponse, "~", paste(MainTerms[!RandTerms], collapse = "+")))

### Check the dispersion model ###
if (is.null(disp)){
	x.disp <- NULL
} else {
	if (!inherits(disp, "formula")) stop("\n Dispersion model must be a one-sided formula of the form \" ~ effect\"")
    if (attr(terms(disp), "response") != 0) stop("\n Dispersion model must be a one-sided formula of the form \" ~ effect\"")
	### Create design matrix for the dispersion model ###
	DispModel <- model.frame(disp, data)
	x.disp <- model.matrix(attr(DispModel, "terms"), data = DispModel)
	row.names(x.disp) <- NULL
}
### Formulae checked ###

### Create design matrix for the fixed effects ###
mf <- match.call(expand.dots = FALSE)
m <- match(c("data", "weights", "offset"), names(mf), 0)
mf <- mf[c(1, m)]
mf$formula <- fixed
mf$drop.unused.levels <- TRUE
mf[[1]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")
Y <- model.response(mf)
X <- if (!is.empty.model(mt)) model.matrix(attr(mf, "terms"), data = mf) else matrix(, NROW(Y), 0)
weights <- as.vector(model.weights(mf))
if (!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
offset <- as.vector(model.offset(mf))
if (!is.null(weights) && any(weights < 0)) stop("negative weights not allowed")
if (!is.null(offset)) {
    if (length(offset) == 1)
        offset <- rep(offset, NROW(Y))
    else if (length(offset) != NROW(Y))
        stop(gettextf("number of offsets is %d should equal %d (number of observations)", length(offset), NROW(Y)), domain = NA)
}
row.names(X) <- NULL
if(is.factor(Y)){
	if(family $ family == "binomial"){
		FactorY <- names(table(Y))
		if(length(FactorY) > 2) warning("More than 2 factors in Binomial response is ambiguous and the last category is considered as success")
		Y <- as.numeric(Y == FactorY[length(FactorY)])
	} else stop(paste("response must be numeric for", family$family, "family."))
}

### Create z matrix ###
RandTerms <- MainTerms[RandTerms]
Z <- NULL
for(i in 1:NrRef) { 
	RanTerm <- unlist(strsplit(RandTerms[i],split="|",fixed=TRUE))
	if(length(RanTerm) == 1) {
		stop(gettextf("Random term %d in mainmodel is not meaningful.", i), domain = NA) 
	} else if (length(RanTerm) == 2) {
		Clust <- get_all_vars(as.formula(paste("~", RanTerm[2])), data = data)
		if (NCOL(Clust) > 1) stop(gettextf("Random term %d in mainmodel contains multiple cluster (grouping) variable.", i), domain = NA) 
		if (NROW(Clust) != NROW(X)) stop(gettextf("Remove all NA before input to the hglm function. Length of cluster/grouping variable in random term %d contradicts.", i), domain = NA)
		Clust <- factor(as.vector(Clust[,1]))
		Col <- as.numeric(unclass(Clust))
		RandLevel <- attributes(Clust)$levels
		if (i == 1) nRand <- length(RandLevel)
		RandCvtmf <- model.frame(as.formula(paste("~", RanTerm[1])), data = data,drop.unused.levels=TRUE)
		RandCvt <- model.matrix(attr(RandCvtmf, "terms"), data = RandCvtmf) # bug fixed by maa 150828
		CheckCatCov <- attributes(RandCvt)$contrasts
		if (!is.null(CheckCatCov)) stop(paste("Categorical covariate",names(CheckCatCov), "not allowed in random effects"))
		if ((NROW(RandCvt) > 0) & (NCOL(RandCvt) > 0)) {
        	for (J in 1:NCOL(RandCvt)){
        		if (is.null(Z)) {
        			Z <- sparseMatrix(i = 1:nrow(data), j = Col, x = as.numeric(RandCvt[,J]), dims = c(nrow(data), length(RandLevel)))
        			colnames(Z) <- paste(colnames(RandCvt)[J], "|", RanTerm[2], ":", RandLevel, sep = "")
        		} else {
        			ZJ <- sparseMatrix(i = 1:nrow(data), j = Col, x = as.numeric(RandCvt[,J]), dims = c(nrow(data), length(RandLevel)))
        			colnames(ZJ) <- paste(colnames(RandCvt)[J], "|", RanTerm[2], ":", RandLevel, sep = "")
        			Z <- cBind(Z, ZJ)
        			nRand <- c(nRand, ncol(ZJ))
        		}
        	}
		}
	} else stop(gettextf("Random term %d in mainmodel contain too many |'s.", i), domain = NA)
}

####Check NA in y, X, Z ####
#### added by lrn 2015-03-24
if ( sum( is.na( model.frame(fixed, na.action=NULL, data=data) ) ) >0 )  warning( "NA in response and/or fixed term. Remove all NA before input to the hglm function.", immediate.=TRUE)
if ( nrow(Z) < nrow( model.frame(fixed, na.action=NULL, data=data ) ) ) warning( "NA in random effects term. Remove all NA before input to the hglm function.", immediate.=TRUE)
if (!is.null(disp)) {
	if ( nrow(x.disp) < nrow( model.frame(fixed, na.action=NULL, data=data ) ) ) warning( "NA in terms of the dispersion model. Remove all NA before input to the hglm function.", immediate.=TRUE)
} 
if (is.null(RandC)) {
	val <- hglm.default(X = X, y = Y, Z = Z, family = family, rand.family = rand.family, X.disp = x.disp,
    	                link.disp = link.disp, method = method, conv = conv, maxit = maxit, startval = startval,
        	            weights = weights, fix.disp = fix.disp, offset = offset, RandC = nRand, sparse = sparse, 
						vcovmat = vcovmat, calc.like = calc.like, bigRR = bigRR, verbose = verbose, ...)
} else {
	val <- hglm.default(X = X, y = Y, Z = Z, family = family, rand.family = rand.family, X.disp = x.disp,
			link.disp = link.disp, method = method, conv = conv, maxit = maxit, startval = startval,
			weights = weights, fix.disp = fix.disp, offset = offset, sparse = sparse, 
			vcovmat = vcovmat, calc.like = calc.like, bigRR = bigRR, verbose = verbose, ...)
}
val$call <- Call

return(val)

}

