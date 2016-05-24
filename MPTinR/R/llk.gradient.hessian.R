
.make.llk.function <- function(model, param.names, length.param.names) {
	model.character <- as.character(unlist(model))
	log.model.character <- lapply(model.character, function(x) paste("log(", x, ")", sep =""))
	data.log.model.character <- lapply(seq_along(model.character), function(x, log.model.character) paste("hank.data.", x, " * ", log.model.character[x], " + ", sep = ""), log.model.character = log.model.character)
	tmp.llk.character <- do.call("paste", as.list(data.log.model.character))
	llk.character <- substr(tmp.llk.character, 0, nchar(tmp.llk.character) - 3)
	parse(text = llk.character)
}

.make.llk.gradient <- function(llk.function, param.names, length.param.names) {
	gradient.functions <- vector("list", length.param.names)
	for (param in seq_along(param.names)) {
		gradient.functions[[param]] <- D(llk.function, param.names[param])
	}
	gradient.functions
}

.make.llk.hessian <- function(llk.function, param.names, length.param.names) {
	gradient.functions <- .make.llk.gradient(llk.function, param.names, length.param.names)
	hessian.functions <- vector("list", length.param.names * length.param.names)
	for (llk.funct.outer in seq_along(gradient.functions)){
		for (llk.funct.inner in seq_along(gradient.functions)){
			hessian.functions[[((llk.funct.outer - 1) * length.param.names) + llk.funct.inner]] <- D(gradient.functions[[llk.funct.outer]], param.names[llk.funct.inner])
		}
	}
	#tmp.out <- matrix(hessian.functions, length.param.names)
	#apply(tmp.out, c(1,2), function(x) parse(text = x[[1]]))
	matrix(hessian.functions, length.param.names)
}
