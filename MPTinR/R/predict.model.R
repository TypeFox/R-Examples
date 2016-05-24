   
gen.predictions <- function(parameter.values, model.filename, restrictions.filename = NULL, n.per.item.type = NULL, model.type = c("easy", "eqn", "eqn2"), reparam.ineq = TRUE, check.model = TRUE){
	
	model.predictions <- function(Q, unlist.model, param.names, n.params, tmp.env){
		#tmpllk.env <- new.env()
		for (i in seq_len(n.params))  assign(param.names[i],Q[i], envir = tmp.env)
		vapply(unlist.model, eval, envir = tmp.env, 0)
	}
	
	##################################
	## functions above, magic below ##
	##################################
	
	if (!is.null(dim(parameter.values))) stop("paramter.values needs to be a vector (i.e., is.null(dim(parameter.values)) == TRUE)!")
	
	tree <- .get.mpt.model(model.filename, model.type, model.check = check.model)
	
	orig.params <- NULL
	use.restrictions <- FALSE
	

	if (!is.null(restrictions.filename)) {
		orig.params <- .find.MPT.params(tree)
		new.restrictions <- .check.restrictions(restrictions.filename, tree)
		if (length(new.restrictions) > 0) use.restrictions <- TRUE
		if (!reparam.ineq) {
			res.no.ineq <- new.restrictions
			for (res in 1:length(new.restrictions)) if (new.restrictions[[res]][3] == "<") res.no.ineq[[1]] <- NULL
			if (length(res.no.ineq) == 0) use.restrictions <- FALSE
			else new.restrictions <- res.no.ineq
			}
		if (use.restrictions) tree <- .apply.MPT.restrictions(tree, new.restrictions)
		restrictions <- new.restrictions
	}
	

	# make arguments:
	
	param.names <- .find.MPT.params(tree)
	length.param.names <- length(param.names)
	categories.per.type <- vapply(tree, length, 0)
	
	if (length.param.names != length(parameter.values)) stop(paste("Length of parameter.values does not correspond to number of model parameters (i.e., model needs ", length.param.names, " parameters, parameter values is of length ", length(parameter.values), ").", sep = ""))
	
	if (all(names(parameter.values) %in% param.names) & all(param.names %in% names(parameter.values))) param.names <- names(parameter.values)
	
	tmpllk.env <- new.env()
	
	predictions <- model.predictions(parameter.values, unlist(tree), param.names, length.param.names, tmp.env = tmpllk.env)

	if (!is.null(n.per.item.type)) {
		if (length(tree) != length(n.per.item.type)) stop(paste("Length of n.per.item.type does not correspond to size of model. Model has ", length(n.per.item.type), " item types (or trees), but n.per.item type is only of length ", length(n.per.item.type), ".", sep = ""))
		n <- vector("numeric", length(unlist(tree)))
		counter <- 1
		for (nt in seq_along(tree)) {
			for (lt in seq_along(tree[[nt]])) {
				n[counter] <- n.per.item.type[nt]
				counter <- counter + 1
			}
		}
		predictions <- predictions * n
	}
	
	return(predictions)
	
}

gen.data <- function(parameter.values, samples, model.filename, data = NULL, n.per.item.type = NULL, restrictions.filename = NULL, model.type = c("easy", "eqn", "eqn2"), reparam.ineq = TRUE, check.model = TRUE){
	
	class.model <- class(model.filename)
	if ("connection" %in% class.model) {
		tmp.model <- readLines(model.filename)
		model.filename <- textConnection(tmp.model)
	}

	tree <- .get.mpt.model(model.filename, model.type, model.check = check.model)
	
	if (is.null(data) & is.null(n.per.item.type)) stop("Either data or n.per.item.type needs to be non-null")
	
	if (!is.null(data) & !is.null(n.per.item.type)) stop("Only one of data or n.per.item.type can be non-null")
	
	if (!is.null(data)) {
		if (!is.vector(data)) stop("data needs to be a vector")
		if (sum(sapply(tree, length)) != length(data)) stop(paste("Size of data does not correspond to size of model (i.e., model needs ", sum(sapply(tree, length)), " datapoints, data gives ", length(data), " datapoints).", sep = ""))
	}
	
	orig.params <- NULL
	use.restrictions <- FALSE
	
	# check if restrictions and or model are connections and if so save them
	class.restr <- class(restrictions.filename)
	if (!is.null(restrictions.filename) & ("connection" %in% class.restr)) {
		tmp.restr <- readLines(restrictions.filename)
		restrictions.filename <- textConnection(tmp.restr)
	}

	if (!is.null(restrictions.filename)) {
		orig.params <- .find.MPT.params(tree)
		new.restrictions <- .check.restrictions(restrictions.filename, tree)
		if (length(new.restrictions) > 0) use.restrictions <- TRUE
		if (!reparam.ineq) {
			res.no.ineq <- new.restrictions
			for (res in 1:length(new.restrictions)) if (new.restrictions[[res]][3] == "<") res.no.ineq[[1]] <- NULL
			if (length(res.no.ineq) == 0) use.restrictions <- FALSE
			else new.restrictions <- res.no.ineq
			}
		if (use.restrictions) tree <- .apply.MPT.restrictions(tree, new.restrictions)
		restrictions <- new.restrictions
	}
	
	param.names <- .find.MPT.params(tree)
	length.param.names <- length(param.names)
	categories.per.type <- vapply(tree, length, 0)
	
	if (is.null(n.per.item.type)) {
		n.per.item.type <- vector("numeric", length(tree))
		for (c in seq_len(length(tree))) {
			n.per.item.type[c] <- sum(data[(sum(categories.per.type[seq_len(c-1)]) + 1):(sum(categories.per.type[seq_len(c)]))])
		}
	}
		
	if (length(tree) != length(n.per.item.type)) stop(paste("Length of n.per.item.type does not correspond to size of model. Model has ", length(tree), " item types (or trees), but n.per.item type is only of length ", length(n.per.item.type), ".", sep = ""))
	
	if (!is.null(restrictions.filename) & ("connection" %in% class.restr)) {
		restrictions.filename <- textConnection(tmp.restr)
	}
	if ("connection" %in% class.model) {
		model.filename <- textConnection(tmp.model)
	}
	
	predictions <- gen.predictions(parameter.values = parameter.values, model.filename = model.filename, restrictions.filename = restrictions.filename, n.per.item.type = NULL, model.type = model.type, reparam.ineq = reparam.ineq, check.model = check.model)
	
	data <- rmultinom(samples, n.per.item.type[1], predictions[1:sum(categories.per.type[1])])
	
	if (length(categories.per.type) > 1) {
		n.data <- vector("list", length(categories.per.type))
		n.data[[1]] <- data
		for (tree in seq_along(categories.per.type)) {
			if (tree == 1) next
			n.data[[tree]] <- rmultinom(samples, n.per.item.type[tree], predictions[(sum(categories.per.type[seq_len(tree-1)])+1):sum(categories.per.type[seq_len(tree)])])		
		}
		data <- do.call("rbind", n.data)
	}
	
	t(data)	
}


sample.data <- function(data, samples, model.filename = NULL, categories.per.type = NULL, model.type = c("easy", "eqn", "eqn2"), check.model = TRUE){
		
	if (!is.vector(data)) stop("data needs to be a vector")
	
	if (is.null(model.filename) & is.null(categories.per.type)) stop("Either mode.filename or categories.per.type needs to be non-null")
	
	if (!is.null(model.filename) & !is.null(categories.per.type)) stop("Only one of mode.filename and categories.per.type can be non-null")
	
	if (!is.null(model.filename)) {
		tree <- .get.mpt.model(model.filename, model.type, model.check = check.model)
		if (sum(sapply(tree, length)) != length(data)) stop(paste("Size of data does not correspond to size of model (i.e., model needs ", sum(sapply(tree, length)), " datapoints, data gives ", length(data), " datapoints).", sep = ""))
		categories.per.type <- vapply(tree, length, 0)
	}
		
	if (sum(categories.per.type) != length(data)) stop("sum(categories.per.type) needs to be equal to length(data).")
	
	n.per.item.type <- vector("numeric", length(categories.per.type))
	for (c in seq_len(length(categories.per.type))) {
		n.per.item.type[c] <- sum(data[(sum(categories.per.type[seq_len(c-1)]) + 1):(sum(categories.per.type[seq_len(c)]))])
	}
	
	predictions <- data / rep(n.per.item.type, times = categories.per.type)
	
	data <- rmultinom(samples, n.per.item.type[1], predictions[1:sum(categories.per.type[1])])
	
	if (length(categories.per.type) > 1) {
		n.data <- vector("list", length(categories.per.type))
		n.data[[1]] <- data
		for (tree in seq_along(categories.per.type)) {
			if (tree == 1) next
			n.data[[tree]] <- rmultinom(samples, n.per.item.type[tree], predictions[(sum(categories.per.type[seq_len(tree-1)])+1):sum(categories.per.type[seq_len(tree)])])		
		}
		data <- do.call("rbind", n.data)
	}
	
	t(data)	
}

