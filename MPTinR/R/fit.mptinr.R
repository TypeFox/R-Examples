
fit.mptinr <- function(data, objective, param.names, categories.per.type, gradient = NULL, use.gradient = TRUE, hessian = NULL, use.hessian = FALSE, prediction = NULL, n.optim = 5, fia.df = NULL, ci = 95, starting.values = NULL, lower.bound = 0, upper.bound = 1, output = c("standard", "fia", "full"), fit.aggregated = TRUE, sort.param = TRUE, show.messages = TRUE, use.restrictions = FALSE, orig.params = NULL, restrictions = NULL, multicore = c("none", "individual", "n.optim"), sfInit = FALSE, nCPU = 2, control = list(), numDeriv = TRUE, ...) {
	
	if (multicore[1] != "none" & sfInit) {
		eval(call("require", package = "snowfall", character.only = TRUE))
        #require(snowfall)
		sfInit( parallel=TRUE, cpus=nCPU )
	} else if (multicore[1] != "none") {
	  if (!eval(call("require", package = "snowfall", character.only = TRUE))) stop("multicore needs snowfall")
	}
	
	n.items.per.type <- function(categories.per.type, data) {
		items.per.type <- array(0, dim = dim(data))
		c1 <- 1
		c2 <- 0
		for (type in categories.per.type) {
		c2 <- c2 + type
		items.per.type[,c1:c2] <- rowSums(data[,c1:c2, drop = FALSE])
		c1 <- c2 + 1
		}
		items.per.type
	}
	
	sat_model <- function(categories.per.type, data){
		NNN <- n.items.per.type(categories.per.type, data)
		temp <- data * log(data/NNN)
		temp[data == 0] <- 0
		llk <- sum(temp)
		return(-llk)
	}
	
	optim.tree <- function(data, objective, gradient, use.gradient, hessian, use.hessian, tmp.env, param.names, n.params, n.optim, start.params, lower.bound, upper.bound, control, ...)  {
		wrapper.nlminb <- function(x, data, objective, gradient, use.gradient, hessian, use.hessian, tmp.env, param.names, n.params, n.optim, start.params, lower.bound, upper.bound, control, ...) {
			if (is.null(start.params)) start.params <- c(0.1, 0.9)
            if (is.list(start.params)) {
                lbound <- start.params[[1]]
                if (length(lbound) == 1) lbound <- rep(lbound, n.params)
                ubound <- start.params[[2]]
                if (length(ubound) == 1) ubound <- rep(ubound, n.params)
                start.params <- vector("numeric", n.params)
                for (c in 1:n.params) {
                start.params[c] <- runif(1, lbound[c], ubound[c])
                }                
            }
			if (length(start.params) == 2) start.params <- runif(n.params, start.params[1], start.params[2])
			nlminb(start.params, objective = objective, gradient = if(use.gradient) gradient, hessian = if(use.hessian) hessian, tmp.env = tmp.env, lower = lower.bound, upper = upper.bound, control = control, data = data, param.names = param.names, n.params = n.params, lower.bound = lower.bound, upper.bound = upper.bound,  ...)
		}
		for (d in seq_along(data)) assign(paste("hank.data.", d, sep = ""), data[d], envir = tmp.env)
		if (multicore[1] == "n.optim") {
			out <- snowfall::sfClusterApplyLB(1:n.optim, wrapper.nlminb, data = data, objective = objective, gradient = gradient, use.gradient = use.gradient, hessian = hessian, use.hessian = use.hessian, tmp.env = tmp.env, param.names = param.names, n.params = n.params, start.params = start.params, lower.bound = lower.bound, upper.bound = upper.bound, control = control, ...)
		} else out <- lapply(1:n.optim, wrapper.nlminb, data = data, objective = objective, gradient = gradient, use.gradient = use.gradient, hessian = hessian, use.hessian = use.hessian, tmp.env = tmp.env, param.names = param.names, n.params = n.params, start.params = start.params, lower.bound = lower.bound, upper.bound = upper.bound, control = control, ...)
		return(out)
	}
	
	optim.mpt <- function(data, objective, gradient, use.gradient, hessian, use.hessian, tmp.env, param.names, n.params, n.optim, start.params, lower.bound, upper.bound, control, n.data, ...) {
		
		minim <- vector("list", n.data)
		data.new <- lapply(1:n.data, function(x, data) data[x,], data = data) 
		llks <- array(NA, dim=c(n.data, n.optim))
		
		if (multicore[1] == "individual" & length(data.new) > 1) {
			 optim.runs <- snowfall::sfClusterApplyLB(data.new, optim.tree, objective = objective, gradient = gradient, use.gradient = use.gradient, hessian = hessian, use.hessian = use.hessian, tmp.env = tmp.env, param.names = param.names, n.params = n.params, n.optim = n.optim, start.params = start.params, lower.bound = lower.bound, upper.bound = upper.bound, control = control, ...)
		} else optim.runs <- lapply(data.new, optim.tree, objective = objective, gradient = gradient, use.gradient = use.gradient, hessian = hessian, use.hessian = use.hessian, tmp.env = tmp.env, param.names = param.names, n.params = n.params, n.optim = n.optim, start.params = start.params, lower.bound = lower.bound, upper.bound = upper.bound, control = control, ...)
		
		for (c.outer in 1:n.data) {
			least.llk <- 1e10
			for (c in 1: n.optim) {
				llks[c.outer, c] <- -(optim.runs[[c.outer]][[c]][["objective"]])
				if (optim.runs[[c.outer]][[c]]["objective"] < least.llk) {
					minim[[c.outer]] <- optim.runs[[c.outer]][[c]]
					least.llk <- optim.runs[[c.outer]][[c]][["objective"]]
				}
			}
		}
		return(list(minim = minim, optim.runs = optim.runs, llks = llks))
	}
	
	get.goodness.of.fit <- function(minim, data, dgf, n.params, n.data, categories.per.type) {
		Log.Likelihood <- sapply(minim, function(x) x$objective)
		G.Squared <- sapply(1:n.data, function(x, data, Log.Likelihood) as.numeric(2*(Log.Likelihood[x]-sat_model(categories.per.type, data[x,,drop = FALSE]))), data = data, Log.Likelihood = Log.Likelihood)
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
			ic <- data.frame(FIA, AIC, BIC, FIA.penalty = fia[,1])
		} else ic <- data.frame(AIC, BIC)
		rownames(ic) <- NULL
		ic
	}
	
	get.model.info <- function(hessian.list, n.params, dgf) {
		rank_hessian <- sapply(hessian.list, function(x) tryCatch(qr(solve(x))$rank, error = function(e) NA))
		return(data.frame(rank.fisher = rank_hessian, n.parameters = n.params, n.independent.categories = dgf))
	}
	
	get.parameter.table.multi <- function(minim, param.names, n.params, n.data, use.restrictions, inv.hess.list, ci, orig.params){
		#recover()
		var.params <- sapply(inv.hess.list, function(x) tryCatch(diag(x), error = function(e) rep(NA, n.params)))
    # as var.params needs to be of nrow=length(parameter.names) (which it isn't if the hessian function fails, the following prohibits termination.
    if (is.null(dim(var.params))) var.params <- matrix(var.params, nrow = length(param.names), ncol = n.data)
		rownames(var.params) <- param.names
		suppressWarnings(confidence.interval <- qnorm(1-((100-ci)/2)/100)*sqrt(var.params))
		estimates <- vapply(minim, "[[", i = "par", minim[[1]]$par)
    if (is.null(dim(estimates))) 
      dim(estimates) <- c(1, length(estimates))
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
					if (restrictions[[c]][3] == "=" & sum(grepl("^\\.?[[:alpha:]]", restrictions[[c]][2]))) {
						parameter_table <- rbind(parameter_table, data.frame(param.names = restrictions[[c]][1], parameter_table[parameter_table$param.names == restrictions[[c]][2], 2:4], restricted.parameter = 1))
					} else if (restrictions[[c]][3] == "=") {
						parameter_table <- rbind(parameter_table, data.frame(param.names = restrictions[[c]][1], estimates = as.numeric(restrictions[[c]][2]), lower.conf = NA, upper.conf = NA, restricted.parameter = 1))
					} else if (restrictions[[c]][3] == "<") {
						tmp.vars <- .find.MPT.params(parse(text = restrictions[[c]][2])[1])
						new.param <- prod(parameter_table.indiv.tmp[parameter_table.indiv.tmp$param.names %in% tmp.vars,2])
						var.tmp <- var.params[rownames(var.params) %in% tmp.vars, counter.n]
						length.var.tmp <- length(var.tmp)
						# Bounds of confindence intervals are computed by the formula given in Baldi & Batchelder (2003, JMP) Equation 19.
						var.bound.tmp <- rep(NA,length.var.tmp)
						for (j in 1:length.var.tmp) var.bound.tmp[j] <- 2*var.tmp[j] + sum(2^(length.var.tmp-1)*(var.tmp[-j]))
						suppressWarnings(ineq.ci <- qnorm(1-((100-ci)/2)/100)*sqrt(min(var.bound.tmp)))
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
        # as var.params needs to be of length(parameter.names) (which it isn't if the hessian function fails, the following prohibits termination.
        if (length(var.params) < length(parameter.names)) var.params <- rep(NA, length(parameter.names))
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
				if (restrictions[[c]][3] == "=" & sum(grepl("^\\.?[[:alpha:]]", restrictions[[c]][2]))) {
					parameter_table <- rbind(parameter_table, data.frame(parameter.names = restrictions[[c]][1], parameter_table[parameter_table$parameter.names == restrictions[[c]][2], 2:4], restricted.parameter = restrictions[[c]][2]))
				} else if (restrictions[[c]][3] == "=") {
					parameter_table <- rbind(parameter_table, data.frame(parameter.names = restrictions[[c]][1], estimates = as.numeric(restrictions[[c]][2]), lower.conf = NA, upper.conf = NA, restricted.parameter = restrictions[[c]][2]))
				} else if (restrictions[[c]][3] == "<") {
					tmp.vars <- .find.MPT.params(parse(text = restrictions[[c]][2])[1])
					new.param <- prod(parameter_table.tmp[parameter_table.tmp$parameter.names %in% tmp.vars,2])
					var.tmp <- var.params[names(var.params) %in% tmp.vars]
					length.var.tmp <- length(var.tmp)
					# Bounds of confindence intervals are computed by the formula given in Baldi & Batchelder (2003, JMP) Equation 19.
					var.bound.tmp <- rep(NA,length.var.tmp)
					for (j in 1:length.var.tmp) var.bound.tmp[j] <- 2*var.tmp[j] + sum(2^(length.var.tmp-1)*(var.tmp[-j]))
					suppressWarnings(ineq.ci <- qnorm(1-((100-ci)/2)/100)*sqrt(min(var.bound.tmp)))
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
	
	get.predicted.values <- function (minim, prediction, data, n_items.per.type, param.names, n.params, n.data, tmp.env, ...) {
		predictions <- array(0, dim = dim(data))
		for (c in 1:n.data){
			for (d in seq_along(data[c,])) assign(paste("hank.data.", d, sep = ""), data[c,d], envir = tmpllk.env)
			tree.eval <- do.call(prediction, args  = list(minim[[c]][["par"]], param.names = param.names, n.params = n.params, tmp.env = tmp.env, ... ))
			frequencies <- n_items.per.type[c,]
			predictions[c,] <- tree.eval * frequencies
		}
		return(predictions)
	}
	
	###############################################################################################
	### above functions for MPTinR, below the code that calls them ################################
	###############################################################################################
	
	#recover()
	
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
	
	if (ci != 95 && show.messages) message(paste("Confidence intervals represent ", ci, "% intervals.", sep = ""))
	
	n.params <- length(param.names)
	length.param.names <- length(param.names)
	
	if (sum(categories.per.type) != length(data[1,])) stop(paste("Size of data does not correspond to size of model (i.e., model needs ", sum(categories.per.type), " datapoints, data gives ", length(data[1,]), " datapoints).", sep = ""))
	
	if (!is.null(starting.values)) {
		if (length(starting.values) != 2) {
			n.optim <- 1
			if (length(starting.values) != n.params) stop("length(starting.values) does not match number of parameters.\nUse check.mpt() to find number and order of parameters!")
		}
	}
	
	if ((length(lower.bound) != 1) && (length(lower.bound) != n.params)) stop("length(lower.bound) does not match number of parameters.\nUse check.mpt() to find number and order of parameters!")
	if ((length(upper.bound) != 1) && (length(upper.bound) != n.params)) stop("length(upper.bound) does not match number of parameters.\nUse check.mpt() to find number and order of parameters!")
	
	if (n.optim != 1 & show.messages) message(paste("Presenting the best result out of ", n.optim, " minimization runs.", sep =""))
	
	if (is.null(fia.df)) fia <- NULL
	else {
		fia <- TRUE
		if (multiFit) {
			fia.df.tmp <- fia.df
			fia.df <- fia.df.tmp[[1]]
			fia.agg.tmp <- fia.df.tmp[[2]]
		}
	}
	
	#n_items <- n.items.per.type(categories.per.type, data)
	n_items <- rowSums(data)
	n_items.per.type <- n.items.per.type(categories.per.type, data)
	n.data <- dim(data)[1]
	dgf <- sum(categories.per.type) - length(categories.per.type)
	
	data.smaller.5 <- t(apply(data, 1, function(x) x < 5))
	# if (any(data.smaller.5)) warning(paste("Categories have n < 5! Do NOT trust these CIs. Dataset:", paste((1:n.data)[apply(data.smaller.5, 1, any)], collapse = " "), sep = ""))
	
	# check if there is a gradient or hessian function:
	if (is.null(gradient)) use.gradient <- FALSE
	if (is.null(hessian)) use.hessian <- FALSE
	
	tmpllk.env <- new.env()
	#attach(tmpllk.env)
	t0 <- Sys.time()
	if (show.messages) {
		print(paste("Model fitting begins at ", t0, sep = ""))
		flush.console()
	}
	res.optim <- optim.mpt(data = data, objective = objective, gradient = gradient, use.gradient = use.gradient, hessian = hessian, use.hessian = use.hessian, tmp.env = tmpllk.env, param.names = param.names, n.params = n.params, n.optim = n.optim, start.params = starting.values, lower.bound = lower.bound, upper.bound = upper.bound, control = control, n.data = n.data, ...)
	
	minim <- res.optim$minim
	optim.runs <- res.optim$optim.runs
	llks <- res.optim$llks
	
	not.converged <- vector("numeric", n.data)
	better.approx <- vector("numeric", n.data)
	better.analytic <- vector("numeric", n.data)
	convergence.second <- vector("numeric", n.data)
	counter.fit.anew <- 1
	mapping.fit.anew <- vector("numeric", n.data)
	
	# the following loop checks if the optimization routine converged succesfully, if not, result is added to another list
	for (c.n in 1:n.data) {
		# warning based on counts not needed, changed to nlminb (HS 26-01-2012)
		#if (minim[[c.n]][["counts"]][1] < 10) warning(paste("Number of iterations run by the optimization routine for individual ", c.n, " is low (i.e., < 10) indicating local minima. Try n.optim >= 5.", sep = ""))
        if (is.null(minim[[c.n]])) stop(paste0("Optimization for dataset ", c.n, " did not evaluate successfully with the given starting values.\nUse appropriate starting values."))
		if (minim[[c.n]][["convergence"]] != 0) {
			not.converged[c.n] <- c.n
			mapping.fit.anew[counter.fit.anew] <- c.n
			counter.fit.anew <- counter.fit.anew + 1
		}
        
	}
	
	if (use.gradient & counter.fit.anew > 1) {
		mapping.fit.anew <- mapping.fit.anew[seq_len(counter.fit.anew-1)]
		new.fits <- suppressWarnings(fit.mptinr(data[mapping.fit.anew,,drop = FALSE], objective = objective, param.names = param.names, categories.per.type = categories.per.type, gradient = gradient, use.gradient = FALSE, hessian = hessian, use.hessian = FALSE, prediction = prediction, n.optim = n.optim, fia.df = NULL, ci = ci, starting.values = starting.values, lower.bound = lower.bound, upper.bound = upper.bound, output = "full", sort.param = sort.param, show.messages = FALSE, use.restrictions = use.restrictions, orig.params = orig.params, restrictions = restrictions, multicore = multicore, sfInit = FALSE, nCPU = 2, control = control, numDeriv = numDeriv, ...))
		if (length(mapping.fit.anew) > 1) {
			for (c.n in seq_along(mapping.fit.anew)) {
				optim.runs[[mapping.fit.anew[c.n]]] <- c(optim.runs[[mapping.fit.anew[c.n]]], new.fits[["optim.runs"]][["individual"]][[c.n]])
				if (n.optim > 1) llks[mapping.fit.anew[c.n],] <- vapply(new.fits[["optim.runs"]][["individual"]][[c.n]], "[[", 0, i = "objective")
				if (new.fits[["best.fits"]][["individual"]][[c.n]][["objective"]] < minim[[mapping.fit.anew[c.n]]][["objective"]]) {
					if (new.fits[["best.fits"]][["individual"]][[c.n]][["convergence"]] != 0) convergence.second[mapping.fit.anew[c.n]] <- new.fits[["best.fits"]][["individual"]][[c.n]][["convergence"]]
					minim[[mapping.fit.anew[c.n]]] <- new.fits[["best.fits"]][["individual"]][[c.n]]
					better.approx[mapping.fit.anew[c.n]] <- mapping.fit.anew[c.n]
				} else {
					better.analytic[mapping.fit.anew[c.n]] <- mapping.fit.anew[c.n]
				}
			}
		} else {
			optim.runs[[mapping.fit.anew[1]]] <- c(optim.runs[[mapping.fit.anew[1]]], new.fits[["optim.runs"]][[1]])
			if (n.optim > 1) llks[mapping.fit.anew[1],] <- vapply(new.fits[["optim.runs"]][[1]], "[[", 0, i = "objective")
			if (new.fits[["best.fits"]][[1]][["objective"]] < minim[[mapping.fit.anew[1]]][["objective"]]) {
				if (new.fits[["best.fits"]][[1]][["convergence"]] != 0) convergence.second[mapping.fit.anew[1]] <- new.fits[["best.fits"]][[1]][["convergence"]]
				minim[[mapping.fit.anew[1]]] <- new.fits[["best.fits"]][[1]]
				better.approx[mapping.fit.anew[1]] <- mapping.fit.anew[1]
			} else {
				better.analytic[mapping.fit.anew[1]] <- mapping.fit.anew[1]
			}
		}
	}

	error.codes <- vapply(minim, "[[", 0, i = "convergence")
	if (sum(not.converged != 0)) {
		if (use.gradient == FALSE) {
			warning(paste("Optimization routine for dataset(s) ", not.converged[not.converged != 0], " did not converge succesfully.
  Error code(s): ", sort(unique(vapply(minim, "[[", 0, i = "convergence")))[-1], ". Try use.gradient == TRUE or use output = 'full' for more information.", sep =""))
		} else {
      if (show.messages) {
  			message(paste("Optimization routine for dataset(s) ", paste(not.converged[not.converged != 0], collapse = " "), "
  did not converge succesfully. Tried again with use.gradient == FALSE.", sep =""))
	  		if (sum(better.approx) != 0) message(paste("Optimization for dataset(s) ", paste(better.approx[better.approx != 0], collapse = " "), "
  using numerically estimated gradients produced better results. Using those results.
  Old results saved in output == 'full' [['optim.runs']].", sep =""))
		  	if (sum(better.analytic) != 0) message(paste("Optimization for dataset(s) ", paste(better.analytic[better.analytic != 0], collapse = " "), "
  using numerical estimated gradients did NOT produce better results.
  Keeping original results. Use output = 'full' for more details.", sep =""))
      }
			if (sum(error.codes) != 0) {
				warning(paste("Error code(s) in final results: ", paste(sort(unique(error.codes)), collapse = " "),  ". The following dataset(s) did not converge succesfully in the best fitting optimization run:
", paste(which(error.codes != 0), collapse = " "), sep =""))
			 }
		}
	}
	best.fits <- minim	
	
	if (show.messages) {
		t1 <- Sys.time()
		print(paste("Model fitting stopped at ", t1, sep = ""))
		print(t1-t0)
	}
	
	hessian.list <- vector("list", n.data)
	
	if (is.null(hessian) & show.messages) message("No function for computing Hessian Matrix specified or it failed. Hessian Matrix is estimated numerically. Validity of CIs is questionable.")
	
	for (c in 1:n.data) {
		for (d in seq_along(data[c,])) assign(paste("hank.data.", d, sep = ""), data[c,d], envir = tmpllk.env)
		if (!is.null(hessian)) hessian.list[[c]] <- tryCatch(do.call(hessian, args  = list(minim[[c]][["par"]], data = data[c,], upper.bound = upper.bound, lower.bound = lower.bound, param.names = param.names, n.params = length.param.names, tmp.env = tmpllk.env, ... )), error = function(e) NA)
		else {
            if (numDeriv) {
              if (show.messages) message("Note: CIs are based on the numerically estimated Hessian matrix")
                hessian.list[[c]] <- tryCatch(numDeriv::hessian(func = objective, x = minim[[c]][["par"]], data = data[c,], upper.bound = upper.bound, lower.bound = lower.bound, param.names = param.names, n.params = length.param.names, tmp.env = tmpllk.env, ... ), error = function(e) NA)
            } else hessian.list[[c]] <- NA
        }
	}
    if (numDeriv && all(vapply(hessian.list, function(x) all(is.na(x)), NA))) {
      if (show.messages) message("Note: CIs are based on the numerically estimated Hessian matrix")
        for (c in 1:n.data) {
            for (d in seq_along(data[c,])) assign(paste("hank.data.", d, sep = ""), data[c,d], envir = tmpllk.env)
            hessian.list[[c]] <- tryCatch(numDeriv::hessian(func = objective, x = minim[[c]][["par"]], data = data[c,], upper.bound = upper.bound, lower.bound = lower.bound, param.names = param.names, n.params = length.param.names, tmp.env = tmpllk.env, ... ), error = function(e) NA)
        }
    }
    
	inv.hess.list <- lapply(hessian.list, function(x) tryCatch(solve(x), error = function(e) NA))
	
	goodness.of.fit <- get.goodness.of.fit(minim, data, dgf, n.params, n.data, categories.per.type = categories.per.type)
	if (is.null(fia)) information.criteria <- get.information.criteria(minim, goodness.of.fit$G.Squared, n.params, n_items)
	else information.criteria <- get.information.criteria(minim, goodness.of.fit$G.Squared, n.params, n_items, fia.df)
	
	model.info <- get.model.info(hessian.list, n.params, dgf)
	if (n.optim > 1) summary.llks <- t(apply(llks, 1, summary))
	
	#recover()
	
	if (multiFit) {
		fia.agg <- NULL
		summed.goodness.of.fit <- data.frame(t(apply(goodness.of.fit, 2, sum)))
		summed.goodness.of.fit[1,4] <- pchisq(summed.goodness.of.fit[1,2], summed.goodness.of.fit[1,3], lower.tail = FALSE)
		goodness.of.fit <- list(individual = goodness.of.fit, sum = summed.goodness.of.fit)
		parameters <- get.parameter.table.multi(minim, param.names, n.params, n.data, use.restrictions, inv.hess.list, ci, orig.params)
    #browser()
    # added correct calculations of BIC sum and FIA sum (March 2015)
    BIC_sum <- summed.goodness.of.fit$G.Squared + (n.params*nrow(data))*log(sum(n_items))
    AIC_sum <- sum(information.criteria$AIC)
  	information.criteria <- list(individual = information.criteria, sum = data.frame(AIC = AIC_sum, BIC = BIC_sum))
    if (!is.null(fia)) {
      #FIA_penalty_sum <- sum(fia.df.tmp[[1]]$lnInt) + sum(fia.df.tmp[[1]]$lnconst) + (n.params*nrow(data))/2*log(sum(n_items)/2/pi)
      FIA_sum <- sum(information.criteria$individual$FIA)
      FIA_penalty_sum <- sum(information.criteria$individual$FIA.penalty)
      information.criteria$sum <- data.frame(FIA = FIA_sum, information.criteria$sum, FIA.penalty = FIA_penalty_sum)
    }
		model.info <- list(individual = model.info)
		if (fit.aggregated) {
			data.pooled <- apply(data,2,sum)
			data.pooled <- matrix(data.pooled, 1, length(data.pooled))
			res.optim.pooled <- optim.mpt(data = data.pooled, objective = objective, gradient = gradient, use.gradient = use.gradient, hessian = hessian, use.hessian = use.hessian, tmp.env = tmpllk.env, param.names = param.names, n.params = n.params, n.optim = n.optim, start.params = starting.values, lower.bound = lower.bound, upper.bound = upper.bound, control = control, n.data = 1, ...)
			
			for (d in seq_along(data.pooled)) assign(paste("hank.data.", d, sep = ""), data.pooled[d], envir = tmpllk.env)
			if (!is.null(hessian)) hessian.pooled <- tryCatch(do.call(hessian, args  = list(res.optim.pooled[["minim"]][[1]][["par"]], upper.bound = upper.bound, lower.bound = lower.bound, data = data.pooled, param.names = param.names, n.params = length.param.names, tmp.env = tmpllk.env, ... )), error = function(e) NA)
			else if (numDeriv) hessian.pooled <- tryCatch(numDeriv::hessian(func = objective, x = res.optim.pooled[["minim"]][[1]][["par"]], upper.bound = upper.bound, lower.bound = lower.bound, data = data.pooled, param.names = param.names, n.params = length.param.names, tmp.env = tmpllk.env, ... ), error = function(e) NA)
            else hessian.pooled <- NA
			if (numDeriv && all(is.na(hessian.pooled))) hessian.pooled <- tryCatch(numDeriv::hessian(func = objective, x = res.optim.pooled[["minim"]][[1]][["par"]], upper.bound = upper.bound, lower.bound = lower.bound, data = data.pooled, param.names = param.names, n.params = length.param.names, tmp.env = tmpllk.env, ... ), error = function(e) NA)
			inv.hessian <- tryCatch(solve(hessian.pooled), error = function(e) NA)
			if (!is.null(fia)) fia.agg <- fia.agg.tmp
			goodness.of.fit <- c(goodness.of.fit, aggregated = list(get.goodness.of.fit(res.optim.pooled[["minim"]], data.pooled, dgf, n.params, 1, categories.per.type = categories.per.type)))
			information.criteria <- c(information.criteria, aggregated = list(get.information.criteria(res.optim.pooled$minim, goodness.of.fit[["aggregated"]][["G.Squared"]], n.params, sum(n_items), fia.agg)))
			model.info <- c(model.info, aggregated = list(get.model.info(list(hessian.pooled), n.params, dgf)))
			parameters <- c(parameters, aggregated = list(get.parameter.table.single(res.optim.pooled[["minim"]][[1]], param.names, n.params, use.restrictions, inv.hessian, ci, sort.param = sort.param, orig.params)))
			if (n.optim > 1) summary.llks <- list(individual = summary.llks, aggregated = summary(res.optim.pooled[["llks"]][[1]]))
			if (output[1] == "full") {
				optim.runs <- c(individual = list(optim.runs), aggregated = res.optim.pooled$optim.runs)
				best.fits <- c(individual = list(minim), aggregated = res.optim.pooled$minim)
				hessian.list <- c(individual = hessian.list, aggregated = list(hessian.pooled))
			}
		} else {
			if (output[1] == "full") {
				optim.runs <- c(individual = list(optim.runs))
				best.fits <- c(individual = list(minim))
				hessian.list <- c(individual = hessian.list)
			}
		}
	} else {
		parameters <- get.parameter.table.single(minim[[1]], param.names, n.params, use.restrictions, inv.hess.list[[1]], ci, sort.param = sort.param, orig.params)
	}
	if (!is.null(prediction)) {
		if (multiFit) {
			predictions <- list(individual = get.predicted.values(minim = minim, prediction = prediction, data = data, tmp.env = tmpllk.env, n_items.per.type = n_items.per.type, param.names = param.names, n.params = n.params, n.data = n.data, ...))
			if (fit.aggregated) {
				n_items.per.type.agg <- n.items.per.type(categories.per.type, data.pooled)
				predictions = c(predictions, aggregated = list(get.predicted.values(minim = res.optim.pooled$minim, prediction = prediction, data = data.pooled, tmp.env = tmpllk.env, n_items.per.type = n_items.per.type.agg, param.names = param.names, n.params = n.params, n.data = 1, ...)))
			}
		} else predictions <- get.predicted.values(minim = minim, prediction = prediction, data = data, tmp.env = tmpllk.env, n_items.per.type = n_items.per.type, param.names = param.names, n.params = n.params, n.data = n.data, ...)
	} else predictions <- list()
	if (multiFit) {
		if (fit.aggregated) data.out <- list(individual = data, aggregated = data.pooled)
		else data.out <- list(individual = data)
	} else data.out <- data
	data <- list(observed = data.out, predicted = predictions)
	
	if (multiFit) if (!is.null(fia.agg)) fia.df <- list(individual = fia.df, aggregated = fia.agg)
	
	outlist <- list(goodness.of.fit = goodness.of.fit, information.criteria = information.criteria, model.info = model.info, parameters = parameters, data = data)
	#if (!is.null(fia)) outlist <- c(outlist, FIA = fia)
	
	if (n.optim > 1) outlist <- c(outlist, fitting.runs = list(summary.llks))
	if (output[1] == "fia" | (output[1] == "full" & !is.null(fia))) outlist <- c(outlist, FIA = list(fia.df))
	if (output[1] == "full") outlist <- c(outlist, optim.runs = list(optim.runs), best.fits = list(best.fits), hessian = list(hessian.list))
	
	if (multicore[1] != "none" & sfInit) snowfall::sfStop()
	return(outlist)
}


