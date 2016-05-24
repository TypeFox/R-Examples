

get.mpt.fia <- function(data, model.filename, restrictions.filename = NULL, Sample = 200000, model.type = c("easy", "eqn", "eqn2"), round.digit = 6, multicore = FALSE, split = NULL, mConst = NULL){
	
	if(is.vector(data)) {
		data <- array(data, dim = c(1, length(data)))
		multiFit <- FALSE
	} else
		if(is.matrix(data) | is.data.frame(data)) {
			if (is.data.frame(data)) data <- as.matrix(data)
			multiFit <- TRUE
	} else stop("data is neither vector, nor matrix, nor data.frame!")
	
	# is model a connection and needs to be reused?
	class.model <- class(model.filename)
	if ("connection" %in% class.model) {
		tmp.model <- readLines(model.filename)
		model.filename <- textConnection(tmp.model)
	}
	
	
	model <- .get.mpt.model(model.filename, model.type)
	n.data <- dim(data)[1]
	
	if(!is.null(restrictions.filename)) {
		restrictions <- .check.restrictions(restrictions.filename, model)
	}
	
	
	
	if (sum(sapply(model, length)) != length(data[1,])) stop(paste("Size of data does not correspond to size of model (i.e., model needs ", sum(sapply(model, length)), " datapoints, data gives ", length(data[1,]), " datapoints).", sep = ""))
	
	df.n <- apply(data, 1, .DF.N.get, tree = model)
	
	n_items <- sapply(df.n, function (x) sum(x[[2]]))
	
	if ("connection" %in% class.model) {
		model.filename <- textConnection(tmp.model)
	}
	mpt.string <- make.mpt.cf(model.filename = model.filename, model.type = model.type)
	
	is.category <- grepl("^[[:digit:]]+$", mpt.string)
	
	s <- paste(ifelse(is.category == 0, "p", "C"), collapse = "")
	
	params <- mpt.string[!is.category]
	category <- mpt.string[is.category]
	
	category <- as.numeric(category)
	
	is.p.join <- grepl("^hank\\.join\\.", params)
	
	c.join <- sum(is.p.join)
	
	p.join <- params[is.p.join]	
	p.n.join <- params[!is.p.join]
	
	hank.restrictions <- vector("list", c.join)
	
	ns <- t(sapply(df.n, function (x) x[[2]]))
	
	p.n.join.lev <- sort(unique(p.n.join))
	f.p.n.join <- factor(p.n.join, levels = p.n.join.lev)
	names(p.n.join.lev) <- 1:length(p.n.join.lev)
	
	if (c.join > 0) {
		for (indiv in 1:n.data) {
			for (c.hank in 1:(c.join)) {
				hank.restrictions[[c.hank]][indiv] <- round(sum(ns[indiv,1:c.hank]) / sum(ns[indiv,1:(c.hank+1)]), round.digit)
			}
		}		
		params.join.mat <- matrix(NA, nrow = n.data, ncol = (c.join))
		for (indiv in 1:n.data) {
			for (par.join in 1:length(p.join)) {
				params.join.mat[indiv,par.join] <- -hank.restrictions[[length(p.join)-(par.join-1)]][indiv]
			}
		}		
		parameters <- vector("list", n.data)
		for (indiv in 1:n.data) {
			parameters[[indiv]] <- c(params.join.mat[indiv,], as.numeric(f.p.n.join))
		}
		
		n.fia <- 1
		i.data <- vector('list', n.data)
		p.fia <- vector("numeric", n.data)
		m.fit <- vector("numeric", n.data)
		i.data[[1]] <- params.join.mat[1,]
		p.fia[1] <- 1
		m.fit[1] <- 1
		if (n.data > 1) {
			for (c in 2:n.data) {
				tmp <- list(params.join.mat[c,])
				if (tmp %in% i.data) p.fia[c] <- which(i.data %in% tmp)
				else {
					n.fia <- n.fia + 1
					i.data[n.fia] <- tmp
					p.fia[c] <- n.fia
					m.fit[n.fia] <- c
				}
			}
		}
	m.fit <- m.fit[1:n.fia]
	} else {
		parameters <- vector("list", n.data)
		for (indiv in 1:n.data) {
			parameters[[indiv]] <- as.numeric(f.p.n.join)
		}
		m.fit <- 1
		p.fia <- rep(1,n.data)
	}	
	if(!is.null(restrictions.filename)) {
		ineq <- vector('list', length(restrictions))
		n.ineq <- 1
		for (restr in 1:length(restrictions)) {
			if (restrictions[[restr]][3] == "=") {
				if (grepl("^[[:digit:]]", restrictions[[restr]][2])) {
					for (indiv in 1:n.data) {
						parameters[[indiv]][params == restrictions[[restr]][1]] <- -as.numeric(restrictions[[restr]][2])
					}
				} else {
					for (indiv in 1:n.data) {
						parameters[[indiv]][params == restrictions[[restr]][1]] <- (parameters[[indiv]][params == restrictions[[restr]][2]])[1]
					}
				}
			} else {
				if (restrictions[[restr]][3] == "<") {
					ineq[[n.ineq]] <- matrix(c((parameters[[indiv]][params == restrictions[[restr]][1]])[1], (parameters[[indiv]][params == restrictions[[restr]][4]])[1]), 1,2)
					n.ineq <- n.ineq + 1
				}
			}
		}
		ineq <- do.call("rbind", ineq)
	} else ineq <- NULL
	
	n.fit <- length(m.fit)
	
	fia.result <- vector('list', n.fit)
	
	for (counter in 1:n.fit) {
		fia.result[[counter]] <- bmpt.fia(s, parameters[[m.fit[counter]]], category, n_items[m.fit[counter]], ineq, Sample, multicore = multicore, split = split, mConst = mConst)
	}
	
	n.params <- length(unique(parameters[[1]][parameters[[1]] > 0]))
	
	res <- vector('list', n.data)
	
	for (c in 1:n.data) {
		res[[c]] <- fia.result[[p.fia[c]]]
		res[[c]][["CFIA"]] <- res[[c]][["lnInt"]] + res[[c]][["lnconst"]]+n.params/2*log(n_items[c]/2/pi)
		res[[c]][["CI.l"]] <- NA
		res[[c]][["CI.u"]] <- NA
	}
	
	as.data.frame(do.call('rbind', res))
	
}


