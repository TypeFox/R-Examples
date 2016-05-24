cvdglars <- function(formula,family=c("binomial","poisson"),data,subset,contrast=NULL,control=list()){
	this.call <- match.call()
	if (missing(data))	data <- environment(formula)
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data", "subset"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	if(attr(mt,"intercept") == 0) stop("Models without intercept are not allowed in this version of the package")
	y <- model.response(mf, "any")
	X <- if (!is.empty.model(mt)) model.matrix(mt,mf,contrasts)
	else stop("Model matrix is empty")
	X <- X[,-1,drop=FALSE]
	fit <- cvdglars.fit(X=X,y=y,family=family,control=control)
	fit$call <- this.call
	fit$formula_cv <- update(formula, as.formula(paste(" ~ ", paste(fit$var_cv,collapse = " + "))))
	fit
}