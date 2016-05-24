
check.mpt <- function(model.filename, restrictions.filename = NULL, model.type = c("easy", "eqn")) {
	if (grepl("\\.eqn$", model.filename) || grepl("\\.EQN$", model.filename)) model.type <- "eqn"
	if (model.type[1] == "eqn") {
		model <- .read.EQN.model(model.filename)
	} else model <- .read.MPT.model(model.filename)
	prob.tree.check <- .check.MPT.probabilities(model)
	if(isTRUE(all(prob.tree.check==1))) {
		prob.corr <- TRUE
	} else {
		prob.corr <- FALSE
		warning(paste("Model not constructed well: Branch probabilities of tree(s) ", paste(which(prob.tree.check!=1), collapse= ", "), " do not sum to 1!", sep = ""))
	}
	orig.params <- .find.MPT.params(model)
	l.orig.params <- length(orig.params)
	n.trees.orig <- length(model)
	orig.model <- list(n.params = l.orig.params, parameters = orig.params)
	parameters <- orig.model
	n.model.categories <- sum(vapply(model, length, 0))
	n.model.df <- n.model.categories - n.trees.orig
	if(!is.null(restrictions.filename)) {
		restrictions <- .read.MPT.restrictions(restrictions.filename)
		model <- .apply.MPT.restrictions(model, restrictions)
		restr.params <- .find.MPT.params(model)
		l.restr.params <- length(restr.params)
		restr.model <- list(n.params = l.restr.params, parameters = restr.params)
		parameters <- list(orig.model = orig.model, restr.model = restr.model)
		message("Inequality restricted parameters are (temporarily) exchanged with dummy parameters called hankX.")
	}
	c(probabilities.eq.1 = prob.corr, n.trees = n.trees.orig, n.model.categories = n.model.categories, n.independent.categories = n.model.df, parameters)
	
}
