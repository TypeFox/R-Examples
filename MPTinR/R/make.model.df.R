
.make.model.df <- function(model) {
	#require(stringr)
	oneLineDeparse <- function(expr){
			paste(deparse(expr), collapse="")
		}
	
	n.trees <- length(model)
	l.trees <- sapply(model, length)
	l.trees <- c(0, l.trees)
	
	fin.model <- vector("list", n.trees)
	
	for (tree in 1:n.trees) {
		utree <- unlist(model[[tree]])
		tree.df.unordered <- do.call("rbind",lapply(1:length(utree), function(t) data.frame(category = t, branches = oneLineDeparse(utree[[t]]), stringsAsFactors = FALSE)))
		
		tree.list <- vector("list", dim(tree.df.unordered)[1])
		for (c1 in 1:length(tree.list)) {
			category <- tree.df.unordered[c1,"category"]
			branch <- strsplit(tree.df.unordered[c1,"branches"], "\\+")
			branch <- gsub(" ", "", branch[[1]])
			tree.list[[c1]] <- data.frame(tree = tree, category = category, branches = branch, stringsAsFactors = FALSE)
		}
		tree.df <- do.call("rbind", tree.list)
		fin.model[[tree]] <- tree.df[rev(order(tree.df[["branches"]])),]
	}
	n.categories <- c(0,sapply(fin.model, function(x) max(x[["category"]])))
	n.cat.cumsum <- cumsum(n.categories)
	
	model.df <- do.call("rbind", fin.model)
	
	model.df[["category"]] <- model.df[,"category"] + n.cat.cumsum[model.df[,"tree"]]
	
	rownames(model.df) <- 1:dim(model.df)[1]
	model.df
	
}
