
make.mpt.cf <- function(model.filename, restrictions.filename = NULL, model.type = c("easy", "eqn"), treewise = FALSE){
	  
	bin.objects <- function(branch) {
		objects <- strsplit(branch, "\\*")[[1]]
		!(grepl("[()]", objects))
	}
	oneLineDeparse <- function(expr){
	  paste(deparse(expr), collapse="")
	}
	model <- .get.mpt.model(model.filename, model.type)
	
  if (!is.null(restrictions.filename)) {
	  new.restrictions <- .check.restrictions(restrictions.filename, model)
	  if (length(new.restrictions) > 0) use.restrictions <- TRUE
    res.no.ineq <- new.restrictions
    for (res in 1:length(new.restrictions)) if (new.restrictions[[res]][3] == "<") res.no.ineq[[1]] <- NULL
    if (length(res.no.ineq) == 0) use.restrictions <- FALSE
    else new.restrictions <- res.no.ineq
	  if (use.restrictions) model <- .apply.MPT.restrictions(model, new.restrictions)
	}

	if (treewise) return(lapply(model, function(x) make.mpt.cf(textConnection(vapply(x, oneLineDeparse, "")))))
  
	model.df <- .make.model.df(model)
	
	#recover()
	
	model.df.tmp <- model.df
	c.join <- 1
	
	while (length(unique(model.df.tmp[,"tree"])) > 1) {
		model.df.tmp[model.df.tmp$tree == 1, "branches"] <- paste("hank.join.", c.join, "*", model.df.tmp[model.df.tmp$tree == 1, "branches"], sep = "")
		model.df.tmp[model.df.tmp$tree == 2, "branches"] <- paste("(1-hank.join.", c.join, ")*", model.df.tmp[model.df.tmp$tree == 2, "branches"], sep = "")
		model.df.tmp[model.df.tmp$tree == 2, "tree"] <- rep(1, length(model.df.tmp[model.df.tmp$tree == 2, "tree"]))
		model.df.tmp[model.df.tmp$tree > 2, "tree"] <- model.df.tmp[model.df.tmp$tree > 2, "tree"] -1
		c.join <- c.join + 1
	}
	tree.ordered <- model.df.tmp
	tree.list <- lapply(1:dim(tree.ordered)[1], function(x) list(category = tree.ordered[x,"category"], branch = tree.ordered[x,"branches"], objects = strsplit(tree.ordered[x,"branches"], "\\*")[[1]], params = .find.MPT.params(tree.ordered[x,"branches"]), binary = bin.objects(tree.ordered[x,"branches"])))
	tmp.tree <- tree.list
	#browser()

	mpt.string <- c(tmp.tree[[1]][["objects"]], tmp.tree[[1]][["category"]])
	for (counter1 in 2:length(tmp.tree)) {
		if (length(tmp.tree[[counter1]][["binary"]]) == length(tmp.tree[[counter1-1]][["binary"]]) & tmp.tree[[counter1-1]][["binary"]][length(tmp.tree[[counter1-1]][["binary"]])] == TRUE & tmp.tree[[counter1]][["binary"]][length(tmp.tree[[counter1]][["binary"]])] == FALSE) {
			mpt.string <- c(mpt.string, tmp.tree[[counter1]][["category"]])
		} else {
		if (length(tmp.tree[[counter1]][["binary"]]) == length(tmp.tree[[counter1-1]][["binary"]]) & tmp.tree[[counter1-1]][["binary"]][length(tmp.tree[[counter1-1]][["binary"]])] == FALSE & tmp.tree[[counter1]][["binary"]][length(tmp.tree[[counter1]][["binary"]])] == TRUE) {
			change <- min(which((tmp.tree[[counter1]][["binary"]] == tmp.tree[[counter1-1]][["binary"]]) == FALSE))+1
			tmp.objects <- tmp.tree[[counter1]][["objects"]][change:(length(tmp.tree[[counter1]][["binary"]]))]
			mpt.string <- c(mpt.string, tmp.objects[tmp.tree[[counter1]][["binary"]][change:length(tmp.tree[[counter1]][["binary"]])]], tmp.tree[[counter1]][["category"]])
		} else {
		if (length(tmp.tree[[counter1]][["binary"]]) > length(tmp.tree[[counter1-1]][["binary"]])) {
			change <- min(which((tmp.tree[[counter1]][["binary"]] == tmp.tree[[counter1-1]][["binary"]][1:length(tmp.tree[[counter1]][["binary"]])]) == FALSE))+1
			if (change < (length(tmp.tree[[counter1-1]][["binary"]]))) {
				tmp.param <- tmp.tree[[counter1]][["objects"]][change:(length(tmp.tree[[counter1]][["binary"]]))]
			} else {
				tmp.new <- tmp.tree[[counter1]][["objects"]][(length(tmp.tree[[counter1-1]][["binary"]])):length(tmp.tree[[counter1]][["binary"]])]
				tmp.param <- tmp.new[tmp.tree[[counter1]][["binary"]][(length(tmp.tree[[counter1-1]][["binary"]])):length(tmp.tree[[counter1]][["binary"]])]]
			}
			mpt.string <- c(mpt.string, tmp.param, tmp.tree[[counter1]][["category"]])
		} else {
		if (length(tmp.tree[[counter1]][["binary"]]) < length(tmp.tree[[counter1-1]][["binary"]])) {
			change <- min(which((tmp.tree[[counter1]][["binary"]] == tmp.tree[[counter1-1]][["binary"]][1:length(tmp.tree[[counter1]][["binary"]])]) == FALSE))+1
			if (change <= length(tmp.tree[[counter1]][["binary"]])) {
			tmp.objects <- tmp.tree[[counter1]][["objects"]][change:(length(tmp.tree[[counter1]][["binary"]]))]
			} else tmp.objects <- NULL
			mpt.string <- c(mpt.string, tmp.objects[tmp.tree[[counter1]][["binary"]][change:length(tmp.tree[[counter1]][["binary"]])]], tmp.tree[[counter1]][["category"]])
		}
		}
		}
		}
	
	}
	mpt.string
}



