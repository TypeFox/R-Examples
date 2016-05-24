
   
fit.mpt.old <- function(data, model.filename, restrictions.filename = NULL, n.optim = 5, fia = NULL, ci = 95, starting.values = NULL, output = c("standard", "fia", "full"), reparam.ineq = TRUE, sort.param = TRUE, model.type = c("easy", "eqn", "eqn2"),  multicore = c("none", "individual", "n.optim"), sfInit = FALSE, nCPU = 2){
	
	if (multicore[1] != "none" & sfInit) {
		eval(call("require", package = "snowfall", character.only = TRUE))
		sfInit( parallel=TRUE, cpus=nCPU )
	} else if (multicore[1] != "none") {
    if (!eval(call("require", package = "snowfall", character.only = TRUE))) stop("multicore needs snowfall")
	}
	
	llk.tree <- function(Q, unlist.tree, data, param.names, length.param.names){
		Q[Q > 1] <- 1
		Q[Q < 0] <- 0
		#tmpllk.env <- new.env()
		for (i in 1:length.param.names)  assign(param.names[i],Q[i], envir = tmpllk.env)
		
		tree.eval <- vapply(unlist.tree, eval, envir = tmpllk.env, 0)
		if (any(tree.eval < 0)) stop(paste("Model not constructed well. Branch (i.e., line) ", which(tree.eval < 0), " produces probabilities < 0!", sep = ""))
		llk <- data * log(tree.eval)
		llk[data == 0] <- 0
		llk <- sum(llk)
	if (is.na(llk)) llk <- -1e10
	if (llk == -Inf) llk <- -1e10
	return(-llk)
	}


	sat_model <- function(tree, data){
		temp.branch <- sapply(tree,length)
		NNN <- rep(.DF.N.get(data,tree)[[2]], temp.branch)
		temp <- data * log(data/NNN)
		temp[data == 0] <- 0
		llk <- sum(temp)
		return(-llk)
	}
	
	optim.tree <- function(data, tree, llk.tree, param.names, n.params, n.optim, method = "L-BFGS-B", start.params)  {
		mpt.optim <- function(x, start.params, llk.tree, tree, data, param.names, n.params, method) {
			if (is.null(start.params)) start.params <- c(0.05, 0.95)
			if (length(start.params) == 2) start.params <- runif(n.params, start.params[1], start.params[2])
			optim(start.params, llk.tree, unlist.tree = tree, data = data, param.names = param.names, length.param.names = n.params, method = method, lower = 0, upper = 1, hessian = TRUE)
		}
		if (multicore[1] == "n.optim") {
			out <- snowfall::sfLapply(1:n.optim, mpt.optim, start.params = start.params, llk.tree = llk.tree, tree = tree, data = data, param.names = param.names, n.params = n.params, method = method)
		} else out <- lapply(1:n.optim, mpt.optim, start.params = start.params, llk.tree = llk.tree, tree = tree, data = data, param.names = param.names, n.params = n.params, method = method)
		return(out)
	}
	
	optim.mpt <- function(data, n.data, tree, llk.tree, param.names, n.params, n.optim, start.params) {
		
		minim <- vector("list", n.data)
		data.new <- lapply(1:n.data, function(x, data) data[x,], data = data) 
		llks <- array(NA, dim=c(n.data, n.optim))
		
		if (multicore[1] == "individual") {
			 optim.runs <- snowfall::sfLapply(data.new, optim.tree, tree = unlist(tree), llk.tree = llk.tree, param.names = param.names, n.params = n.params, n.optim = n.optim, start.params = start.params)
		} else optim.runs <- lapply(data.new, optim.tree, tree = unlist(tree), llk.tree = llk.tree, param.names = param.names, n.params = n.params, n.optim = n.optim, start.params = start.params)
		
		for (c.outer in 1:n.data) {
			least.llk <- 1e10
			for (c in 1: n.optim) {
				llks[c.outer, c] <- -(optim.runs[[c.outer]][[c]][["value"]])
				if (optim.runs[[c.outer]][[c]]["value"] < least.llk) {
					minim[[c.outer]] <- optim.runs[[c.outer]][[c]]
					least.llk <- optim.runs[[c.outer]][[c]][["value"]]
				}
			}
		}
		return(list(minim = minim, optim.runs = optim.runs, llks = llks))
	}
	
	get.goodness.of.fit <- function(minim, tree, data, dgf, n.params, n.data) {
		Log.Likelihood <- sapply(minim, function(x) x$value)
		G.Squared <- sapply(1:n.data, function(x, data, Log.Likelihood) as.numeric(2*(Log.Likelihood[x]-sat_model(tree, data[x,]))), data = data, Log.Likelihood = Log.Likelihood)
		df <- dgf - n.params
		p.value <- pchisq(G.Squared,df,lower.tail=FALSE)
		data.frame(Log.Likelihood = -Log.Likelihood, G.Squared, df, p.value)
	}
	
	get.information.criteria <- function(minim, G.Squared, n.params, n_items, fia = NULL) {
		AIC <- G.Squared + 2*n.params
		BIC <- G.Squared + n.params*log(n_items)
		
		#The following three lines would calculate ICOMP, but it is currently deactivated
		#diag.hess <- sapply(inv.hess.list, function(x) tryCatch(diag(x), error = function(e) rep(0, n.params)))
		#det.hess <- sapply(inv.hess.list, function(x) tryCatch(det(x), error = function(e) NA))
		#ICOMP <- G.Squared + n.params*log(colSums(diag.hess)/n.params)-log(det.hess)
		
		if (!is.null(fia)) {
			FIA <- (G.Squared/2) + fia[,1]
			ic <- data.frame(FIA, AIC, BIC)
		} else ic <- data.frame(AIC, BIC)
		rownames(ic) <- NULL
		ic
	}
	
	get.model.info <- function(minim, n.params, dgf) {
		rank_hessian <- sapply(minim, function(x) qr(x$hessian)$rank)
		return(data.frame(rank.hessian = rank_hessian, n.parameters = n.params, n.independent.categories = dgf))
	}
	
	get.parameter.table.multi <- function(minim, param.names, n.params, n.data, use.restrictions, inv.hess.list, ci, orig.params){
		#recover()
		var.params <- sapply(inv.hess.list, function(x) tryCatch(diag(x), error = function(e) rep(NA, n.params)))
		rownames(var.params) <- param.names
		confidence.interval <- qnorm(1-((100-ci)/2)/100)*sqrt(var.params)
		estimates <- sapply(minim, function(x) x$par)
		upper.conf <- estimates + confidence.interval
		lower.conf <- estimates - confidence.interval
		
		if(!use.restrictions) {
			params <- 1:n.params
			names(params) <- param.names
			tmp.values <- NULL
			for (counter.n in 1:n.data) {
				tmp.values <- c(tmp.values, estimates[, counter.n], lower.conf[, counter.n], upper.conf[, counter.n])
			}
			parameter_array <- array(tmp.values, dim = c(n.params, 3, n.data))
			dimnames(parameter_array) <- list(param.names, c("estimates", "lower.conf", "upper.conf"), paste("dataset:", 1:n.data))
			mean.df = data.frame(estimates = apply(parameter_array, c(1,2), mean)[,1], lower.conf = NA, upper.conf = NA)
			rownames(mean.df) <- param.names
		}
		
		if(use.restrictions) {
			
			#parameter_table.tmp <- data.frame(param.names, estimates, lower.conf, upper.conf, restricted.parameter = "")
			
			used.rows <- param.names %in% orig.params
			params <- 1:length(orig.params)
			parameter.names.all <- param.names[used.rows]
			restricted <- rep("", sum(used.rows))
			for (c in 1:length(restrictions)) {
				parameter.names.all <- c(parameter.names.all, restrictions[[c]][1])
				restricted <- c(restricted, restrictions[[c]][3])
			}
			names(params) <- parameter.names.all
			pnames <- parameter.names.all
			
			tmp.values <- NULL
			for (counter.n in 1:n.data) {
				parameter_table.indiv.tmp <- data.frame(param.names, estimates = estimates[, counter.n], lower.conf =lower.conf[, counter.n], upper.conf = upper.conf[, counter.n], restricted.parameter = 0)
				parameter_table <- parameter_table.indiv.tmp[parameter_table.indiv.tmp$param.names %in% orig.params,]
				for (c in 1:length(restrictions)) {
					if (restrictions[[c]][3] == "=" & sum(grepl("[[:alpha:]]", restrictions[[c]][2]))) {
						parameter_table <- rbind(parameter_table, data.frame(param.names = restrictions[[c]][1], parameter_table[parameter_table$param.names == restrictions[[c]][2], 2:4], restricted.parameter = 1))
					}
					if (restrictions[[c]][3] == "=" & sum(grepl("^[[:digit:]]\\.?[[:digit:]]*", restrictions[[c]][2]))) {
						parameter_table <- rbind(parameter_table, data.frame(param.names = restrictions[[c]][1], estimates = as.numeric(restrictions[[c]][2]), lower.conf = NA, upper.conf = NA, restricted.parameter = 1))
					}
					if (restrictions[[c]][3] == "<") {
						tmp.vars <- .find.MPT.params(parse(text = restrictions[[c]][2])[1])
						new.param <- prod(parameter_table.indiv.tmp[parameter_table.indiv.tmp$param.names %in% tmp.vars,2])
						var.tmp <- var.params[rownames(var.params) %in% tmp.vars, counter.n]
						length.var.tmp <- length(var.tmp)
						# Bounds of confindence intervals are computed by the formula given in Baldi & Batchelder (2003, JMP) Equation 19.
						var.bound.tmp <- rep(NA,length.var.tmp)
						for (j in 1:length.var.tmp) var.bound.tmp[j] <- 2*var.tmp[j] + sum(2^(length.var.tmp-1)*(var.tmp[-j]))
						ineq.ci <- qnorm(1-((100-ci)/2)/100)*sqrt(min(var.bound.tmp))
						parameter_table <- rbind(parameter_table, data.frame(param.names = restrictions[[c]][1], estimates = new.param, lower.conf = new.param - ineq.ci, upper.conf = new.param + ineq.ci, restricted.parameter = 2))
					}
				}
				rownames(parameter_table) <- parameter_table$param.names
				parameter_table <- parameter_table[,-1]
				if (sort.param) parameter_table <- parameter_table[order(names(params)),]
				tmp.values <- c(tmp.values, as.vector(as.matrix(parameter_table)))
			}
			if (sort.param) {
				order.new <- order(pnames)
				pnames <- pnames[order.new]
				restricted <- restricted[order.new]
			}
			parameter_array <- array(tmp.values, dim = c(length(orig.params), 4, n.data))
			dimnames(parameter_array) <- list(pnames, c("estimates", "lower.conf", "upper.conf", "restricted.parameter"), paste("dataset:", 1:n.data))
			mean.df = data.frame(apply(parameter_array, c(1,2), mean)[,1], lower.conf = NA, upper.conf = NA, restricted = restricted)
			rownames(mean.df) <- pnames
			colnames(mean.df) <- c("estimates", "lower.conf", "upper.conf", "restricted.parameter")
		}
		return(list(individual = parameter_array, mean = mean.df))
	}

	get.parameter.table.single <- function(minim, parameter.names, n.params, use.restrictions, inv.hess, ci, sort.param, orig.params){
		var.params <- tryCatch(diag(inv.hess), error = function(e) rep(NA, n.params))
		names(var.params) <- parameter.names
		confidence.interval <- tryCatch(qnorm(1-((100-ci)/2)/100)*sqrt(var.params), error = function(e) rep(NA, n.params))
		
		estimates <- minim$par
		upper.conf <- estimates + confidence.interval
		lower.conf <- estimates - confidence.interval
		param.names <- parameter.names
		if(!use.restrictions) parameter_table <- data.frame(estimates, lower.conf, upper.conf)
		if(use.restrictions) {
			parameter_table.tmp <- data.frame(parameter.names, estimates, lower.conf, upper.conf, restricted.parameter = "")
			
			parameter_table <- parameter_table.tmp[parameter_table.tmp$parameter.names %in% orig.params,]
			
			for (c in 1:length(restrictions)) {
				if (restrictions[[c]][3] == "=" & sum(grepl("[[:alpha:]]", restrictions[[c]][2]))) {
					parameter_table <- rbind(parameter_table, data.frame(parameter.names = restrictions[[c]][1], parameter_table[parameter_table$parameter.names == restrictions[[c]][2], 2:4], restricted.parameter = restrictions[[c]][2]))
				}
				if (restrictions[[c]][3] == "=" & sum(grepl("^[[:digit:]]\\.?[[:digit:]]*", restrictions[[c]][2]))) {
					parameter_table <- rbind(parameter_table, data.frame(parameter.names = restrictions[[c]][1], estimates = as.numeric(restrictions[[c]][2]), lower.conf = NA, upper.conf = NA, restricted.parameter = restrictions[[c]][2]))
				}
				if (restrictions[[c]][3] == "<") {
					tmp.vars <- .find.MPT.params(parse(text = restrictions[[c]][2])[1])
					new.param <- prod(parameter_table.tmp[parameter_table.tmp$parameter.names %in% tmp.vars,2])
					var.tmp <- var.params[names(var.params) %in% tmp.vars]
					length.var.tmp <- length(var.tmp)
					# Bounds of confindence intervals are computed by the formula given in Baldi & Batchelder (2003, JMP) Equation 19.
					var.bound.tmp <- rep(NA,length.var.tmp)
					for (j in 1:length.var.tmp) var.bound.tmp[j] <- 2*var.tmp[j] + sum(2^(length.var.tmp-1)*(var.tmp[-j]))
					ineq.ci <- qnorm(1-((100-ci)/2)/100)*sqrt(min(var.bound.tmp))
					parameter_table <- rbind(parameter_table, data.frame(parameter.names = restrictions[[c]][1], estimates = new.param, lower.conf = new.param - ineq.ci, upper.conf = new.param + ineq.ci, restricted.parameter = "<"))
				}
			}
			param.names <- as.character(parameter_table[,1])
			parameter_table <- parameter_table[,2:5]
		}
		if (sort.param) {
			parameter_table <- parameter_table[order(param.names),]
			param.names <- param.names[order(param.names)]
		}
		rownames(parameter_table) <- param.names
		return(parameter_table)
	}
	
	get.predicted.values <- function (minim, tree, data, df.n, param.names, length.param.names, n.data) {
		predictions <- matrix(NA, n.data, length(data[1,]))
		predict.env <- new.env()
		temp.branch <- sapply(tree,length)
		for (c in 1:n.data){
			for (i in 1:length.param.names)  assign(param.names[i],minim[[c]][["par"]][i], envir = predict.env)
			tree.eval <- sapply(unlist(tree), eval, envir = predict.env)
			frequencies <- rep(df.n[[c]][[2]], temp.branch)
			predictions[c,] <- tree.eval * frequencies
		}
		return(predictions)
	}
	
	###############################################################################################
	### above functions for MPTinR, below the code that calls them ################################
	###############################################################################################
	
	#recover()
	
	tree <- .get.mpt.model(model.filename, model.type)
	if(is.null(data)) stop("Model seems to be constructed well (i.e., all probabilities sum to 1), but data is NULL.")
	
	if(is.vector(data)) {
		data <- array(data, dim = c(1, length(data)))
		multiFit <- FALSE
	} else
		if(dim(data)[1] == 1) {
			if (is.data.frame(data)) data <- as.matrix(data)
			data <- array(data, dim = c(1,length(data)))
			multiFit <- FALSE
		} else 
			if(is.matrix(data) | is.data.frame(data)) {
				if (is.data.frame(data)) data <- as.matrix(data)
				multiFit <- TRUE
			} else stop("data is neither vector, nor matrix, nor data.frame!")
	
	if (sum(sapply(tree, length)) != length(data[1,])) stop(paste("Size of data does not correspond to size of model (i.e., model needs ", sum(sapply(tree, length)), " datapoints, data gives ", length(data[1,]), " datapoints).", sep = ""))
	
	if (ci != 95) message(paste("Confidence intervals represent ", ci, "% intervals.", sep = ""))
	
	orig.params <- NULL
	use.restrictions <- FALSE
	
	if (!is.null(restrictions.filename)) {
		use.restrictions <- TRUE
		restrictions <- .read.MPT.restrictions(restrictions.filename)
		orig.tree <- tree
		orig.params <- .find.MPT.params(tree)
		if (!reparam.ineq) {
			res.no.ineq <- restrictions
			for (res in 1:length(restrictions)) if (restrictions[[res]][3] == "<") res.no.ineq[[1]] <- NULL
			if (length(res.no.ineq) == 0) use.restrictions <- FALSE
			else restrictions <- res.no.ineq
			}
		if (use.restrictions) tree <- .apply.MPT.restrictions(tree, restrictions)
	}
	
	param.names <- .find.MPT.params(tree)
	n.params <- length(param.names)
	
	
	if (!is.null(starting.values)) {
		if (length(starting.values) != 2) {
			n.optim <- 1
			if (length(starting.values) != n.params) stop("length(starting.values) does not match number of parameters.\nUse check.mpt() to find number and order of parameters!")
		}
	}
	
	if (n.optim != 1) message(paste("Presenting the best result out of ", n.optim, " minimization runs.", sep =""))
	
	df.n <- apply(data, 1, .DF.N.get, tree = tree)
	n_items <- sapply(df.n, function (x) sum(x[[2]]))
	dgf <- df.n[[1]][[1]]
	n.data <- dim(data)[1]
	
	data.smaller.5 <- t(apply(data, 1, function(x) x < 5))
	if (any(data.smaller.5)) warning(paste("Categories have n < 5! Do NOT trust these CIs. Dataset:", paste((1:n.data)[apply(data.smaller.5, 1, any)], collapse = " "), sep = ""))
	
	tmpllk.env <- new.env()
	t0 <- Sys.time()
	print(paste("Model fitting begins at ", t0, sep = ""))
	flush.console()
	res.optim <- optim.mpt(data, n.data, tree, llk.tree, param.names, n.params, n.optim, starting.values)
	t1 <- Sys.time()
	print(paste("Model fitting stopped at ", t1, sep = ""))
	print(t1-t0)
	
	minim <- res.optim$minim
	optim.runs <- res.optim$optim.runs
	llks <- res.optim$llks
	
	if (!is.null(fia)) {
		if (multiFit) {
			data.new <- rbind(data, apply(data,2,sum))
			fia.tmp <- get.mpt.fia(data.new, model.filename, restrictions.filename, fia, model.type)
			fia.df <- fia.tmp[-dim(fia.tmp)[1],]
			fia.agg.tmp <- fia.tmp[dim(fia.tmp)[1],]
		} else {
			fia.df <- get.mpt.fia(data, model.filename, restrictions.filename, fia, model.type)
		}
	}
	
	inv.hess.list <- lapply(minim, function(x) tryCatch(solve(x$hessian), error = function(e) NA))
	
	goodness.of.fit <- get.goodness.of.fit(minim,tree, data, dgf, n.params, n.data)
	if (is.null(fia)) information.criteria <- get.information.criteria(minim, goodness.of.fit$G.Squared, n.params, n_items)
	else information.criteria <- get.information.criteria(minim, goodness.of.fit$G.Squared, n.params, n_items, fia.df)
	
	model.info <- get.model.info(minim, n.params, dgf)
	if (n.optim > 1) summary.llks <- t(apply(llks, 1, summary))
	
	#recover()
	
	if (multiFit) {
		fia.agg <- NULL
		data.pooled <- apply(data,2,sum)
		data.pooled <- matrix(data.pooled, 1, length(data.pooled))
		res.optim.pooled <- optim.mpt(data.pooled, 1, tree, llk.tree, param.names, n.params, n.optim, starting.values)
		inv.hessian <- tryCatch(solve(res.optim.pooled[["minim"]][[1]][["hessian"]]), error = function(e) NA)
		if (!is.null(fia)) fia.agg <- fia.agg.tmp
		summed.goodness.of.fit <- data.frame(t(apply(goodness.of.fit, 2, sum)))
		summed.goodness.of.fit[1,4] <- pchisq(summed.goodness.of.fit[1,2], summed.goodness.of.fit[1,3], lower.tail = FALSE)
		goodness.of.fit <- list(individual = goodness.of.fit, sum = summed.goodness.of.fit, aggregated = get.goodness.of.fit(res.optim.pooled[["minim"]], tree, data.pooled, dgf, n.params, 1))
		information.criteria <- list(individual = information.criteria, sum = data.frame(t(apply(information.criteria, 2, sum))), aggregated = get.information.criteria(res.optim.pooled$minim, goodness.of.fit[["aggregated"]][["G.Squared"]], n.params, sum(n_items), fia.agg))
		model.info <- list(individual = model.info, aggregated = get.model.info(res.optim.pooled$minim, n.params, dgf))
		parameters <- c(get.parameter.table.multi(minim, param.names, n.params, n.data, use.restrictions, inv.hess.list, ci, orig.params), aggregated = list(get.parameter.table.single(res.optim.pooled[["minim"]][[1]], param.names, n.params, use.restrictions, inv.hessian, ci, sort.param = sort.param, orig.params)))
		if (n.optim > 1) summary.llks <- list(individual = summary.llks, aggregated = summary(res.optim.pooled[["llks"]][[1]]))
		if (output[1] == "full") optim.runs <- c(individual = list(optim.runs), aggregated = res.optim.pooled$optim.runs)
		for (c.n in 1:n.data) {
			if (minim[[c.n]][["counts"]][1] < 10) warning(paste("Number of iterations run by the optimization routine for individual ", c.n, " is low (i.e., < 10) indicating local minima. Try n.optim >= 5.", sep = ""))
			if (minim[[c.n]][["convergence"]] != 0) warning(paste("Optimization routine for individual ", c.n, " did not converge succesfully. Error code: ", minim[[c.n]][["convergence"]], ". Use output = 'full' for more information.", sep =""))
		}
	} else {
		parameters <- get.parameter.table.single(minim[[1]], param.names, n.params, use.restrictions, inv.hess.list[[1]], ci, sort.param = sort.param, orig.params)
		if (minim[[1]][["counts"]][1] < 10) warning("Number of iterations run by the optimization routine is low (i.e., < 10) indicating local minima. Try n.optim >= 5.")
		if (minim[[1]][["convergence"]] != 0) warning(paste("Optimization routine did not converge succesfully. Error code is, ", minim[[1]][["convergence"]], ". Use output = 'full' for more information.", sep =""))
	}
	predictions <- get.predicted.values(minim, tree, data, df.n, param.names, n.params, n.data)
	data <- list(observed = data, predicted = predictions)
	
	if (multiFit) if (!is.null(fia.agg)) fia.df <- list(individual = fia.df, aggregated = fia.agg)
	
	
	outlist <- list(goodness.of.fit = goodness.of.fit, information.criteria = information.criteria, model.info = model.info, parameters = parameters, data = data)
	#if (!is.null(fia)) outlist <- c(outlist, FIA = fia)
	
	if (n.optim > 1) outlist <- c(outlist, summary.llks = list(summary.llks))
	if (output[1] == "fia" | (output[1] == "full" & !is.null(fia))) outlist <- c(outlist, FIA = list(fia.df))
	if (output[1] == "full") outlist <- c(outlist, optim.runs = list(optim.runs))
	
	if (multicore[1] != "none" & sfInit) snowfall::sfStop()
	return(outlist)
}



