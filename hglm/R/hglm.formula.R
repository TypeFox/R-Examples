`hglm.formula` <-
	function(X = NULL, y = NULL, Z = NULL, family = gaussian(link = identity),
			rand.family = gaussian(link = identity), method = "EQL", conv = 1e-6, maxit = 50, 
			startval = NULL, fixed = NULL, random = NULL, X.disp = NULL, disp = NULL, 
			link.disp = "log", X.rand.disp = NULL, rand.disp = NULL, link.rand.disp = "log", 
			data = NULL, weights = NULL, fix.disp = NULL, offset = NULL, RandC = ncol(Z), 
			sparse = TRUE, vcovmat = FALSE, calc.like = FALSE, bigRR = FALSE, verbose = FALSE, ...) {

Call <- match.call()

### check fixed effects formula ###
if (!inherits(fixed, "formula") || length(fixed) != 3) stop("\n Fixed-effects model must be a formula of the form \"resp ~ pred\"")

### Get GLM family and link ###
### Check random effects ###
if (!inherits(random, "formula")) stop("\n Random part must be a one-sided formula of the form \" ~ effect|Cluster\"")
if (attr(terms(random), "response") != 0) stop("\n Random part must be a one-sided formula of the form \" ~ effect|Cluster\"")
if (all.names(random)[2] != "|") stop("The subjects/clusters in Random must be separated by \" |\"")

### Check the dispersion model ###
if (is.null(disp)) {
    x.disp <- NULL
} else {
	if (!inherits(disp, "formula")) stop("\n Dispersion model must be a one-sided formula of the form \" ~ effect\"")
    if (attr(terms(disp), "response") != 0) stop("\n Dispersion model must be a one-sided formula of the form \" ~ effect\"")
	### Create design matrix for the dispersion model ###
	DispModel <- model.frame(disp, data)
	x.disp <- model.matrix(attr(DispModel, "terms"), data = DispModel)
	row.names(x.disp) <- NULL
}
### random effects part is checked ###

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
if (is.factor(Y)) {
	if (family$family == "binomial") {
		FactorY <- names(table(Y))
		if (length(FactorY) > 2) warning("More than 2 factors in Binomial response is ambiguous and the last category is considered as success")
		Y <- as.numeric(Y == FactorY[length(FactorY)])
	} else stop(paste("response must be numeric for", family$family, "family."))
}

### Create z matrix ###
RanTerm <- unlist(strsplit(attr(terms(random), "term.labels"), split = "|", fixed = TRUE))
if (length(RanTerm) > 2) stop("Currently only one random term is supported for hglm(). Consider using hglm2().")
RanTerm <- gsub(pattern = " ", replacement = "", RanTerm)
if (!is.factor(data[1:2, RanTerm[2]])) {
	if ((length(RanTerm) == 2) & (RanTerm[1] == "1")) {     
		ranf <- paste("~", "as.factor(", RanTerm[2], ")", "-1", sep = "")
	} else {
		ranf <- paste("~", "(", RanTerm[1], ")", ":as.factor(", RanTerm[2], ")", "-1", sep = "")
	}
} else {
	if ((length(RanTerm) == 2) & (RanTerm[1] == "1")) {
		ranf <- paste("~", RanTerm[2], "-1", sep = "")
	} else {
		ranf <- paste("~", "(", RanTerm[1], ")", ":", RanTerm[2], "-1", sep = "")
    }
}
ranf <- as.formula(ranf)
rmf <- model.frame(ranf, data, drop.unused.levels = TRUE) # bug fixed by maa 150828
z <- model.matrix(attr(rmf, "terms"), data = rmf)
row.names(z) <- NULL

####Check NA in y, X, Z ####
#### added by lrn 2015-03-24
if ( sum(is.na( model.frame(fixed, na.action=NULL, data=data))) > 0) warning( "NA in response and/or fixed term. Remove all NA before input to the hglm function.", immediate.=TRUE)
if ( sum(is.na( model.frame(ranf, na.action=NULL, data=data))) > 0) warning( "NA in random effects term. Remove all NA before input to the hglm function.", immediate.=TRUE)
if (!is.null(disp)) {
	if ( nrow(x.disp) < nrow( model.frame(fixed, na.action=NULL, data=data ) ) ) warning( "NA in terms of the dispersion model. Remove all NA before input to the hglm function.", immediate.=TRUE)
} 

val <- hglm.default(X = X, y = Y, Z = z, family = family, rand.family = rand.family, X.disp = x.disp,
                    link.disp = link.disp, method = method, conv = conv, maxit = maxit, startval = startval,
                    weights = weights, fix.disp = fix.disp, offset = offset, sparse = sparse, vcovmat = vcovmat, 
					calc.like = calc.like, bigRR = bigRR, verbose = verbose, ...)
val$call <- Call

return(val)

}

