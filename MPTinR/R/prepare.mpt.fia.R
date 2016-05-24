
prepare.mpt.fia <- function(data, model.filename, restrictions.filename = NULL, outfile = "clipboard", Sample = 200000, model.type = c("easy", "eqn", "eqn2")){
	
	if(is.vector(data)) {
		data <- array(data, dim = c(1, length(data)))
		multiFit <- FALSE
	} else
		if(is.matrix(data) | is.data.frame(data)) {
			if (is.data.frame(data)) data <- as.matrix(data)
			multiFit <- TRUE
	} else stop("data is neither vector, nor matrix, nor data.frame!")
	
	if(!is.null(restrictions.filename)) {
		restrictions <- .read.MPT.restrictions(restrictions.filename)
	}
  class.model <- class(model.filename)
  if ("connection" %in% class.model) {
    tmp.model <- readLines(model.filename)
    model.filename <- textConnection(tmp.model)
    model <- .get.mpt.model(model.filename, model.type)
    model.filename <- textConnection(tmp.model)
  } else model <- .get.mpt.model(model.filename, model.type)
	
	n.data <- dim(data)[1]
	
	if (sum(sapply(model, length)) != length(data[1,])) stop(paste("Size of data does not correspond to size of model (i.e., model needs ", sum(sapply(model, length)), " datapoints, data gives ", length(data[1,]), " datapoints).", sep = ""))
	
	df.n <- apply(data, 1, .DF.N.get, tree = model)
	
	n_items <- sapply(df.n, function (x) sum(x[[2]]))
	
	mpt.string <- make.mpt.cf(model.filename, model.type=model.type)
	
	is.category <- grepl("^[[:digit:]]+$", mpt.string)
	
	s <- paste(ifelse(is.category == 0, "p", "C"), collapse = "")
	
	params <- mpt.string[!is.category]
	category <- mpt.string[is.category]
	
	category <- paste(category, collapse = ",")
	
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
				hank.restrictions[[c.hank]][indiv] <- sum(ns[indiv,1:c.hank]) / sum(ns[indiv,1:(c.hank+1)])
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
		parameters[[indiv]] <- c(as.character(round(params.join.mat[indiv,],4)), as.numeric(f.p.n.join))
	}
	} else {
		parameters <- vector("list", n.data)
		for (indiv in 1:n.data) {
			parameters[[indiv]] <- as.character(as.numeric(f.p.n.join))
		}
	}
	
	ineq <- "["
	if(!is.null(restrictions.filename)) {
		for (restr in 1:length(restrictions)) {
			if (restrictions[[restr]][3] == "=") {
				if (grepl("^[[:digit:]]", restrictions[[restr]][2])) {
					for (indiv in 1:n.data) {
						parameters[[indiv]][params == restrictions[[restr]][1]] <- as.character(-as.numeric(restrictions[[restr]][2]))
					}
				} else {
					for (indiv in 1:n.data) {
						parameters[[indiv]][params == restrictions[[restr]][1]] <- (parameters[[indiv]][params == restrictions[[restr]][2]])[1]
					}
				}
			} else {
				if (restrictions[[restr]][3] == "<") {
					if (nchar(ineq) > 1) ineq <- paste(ineq, ";", sep = "")
					ineq <- paste(ineq, paste((parameters[[indiv]][params == restrictions[[restr]][1]])[1], (parameters[[indiv]][params == restrictions[[restr]][4]])[1], sep = ","), sep = "")
				}
			}
		}
	}
	
	ineq <- paste(ineq, "]", sep = "")
	
	calls <- vector("character", n.data)
	
	for (indiv in 1:n.data) {
		calls[indiv] <- paste("[CFIA,CI,lnInt,CI1,lnconst,CI2] = BMPTFIA(\'", s, "\',[", paste(parameters[[indiv]], collapse=","), "],", ineq, ",[", category, "],", n_items[indiv], ",", Sample, ")", sep ="")
	}
	
	writeLines(calls, outfile)
	
	outlist <- list(s = s, parameters = parameters, param.codes = p.n.join.lev, category = category, ineq0 = ineq, n = n_items, internal = mpt.string)
	
	outlist

}



