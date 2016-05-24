# last modified 2014-10-27 by J. Fox

combineModels <- function (...){
	UseMethod("combineModels")
}

combineModels.semmod <- function(..., warn=TRUE){
	model <- rbind(...)
	removeRedundantPaths(model, warn=warn)
	model
}

update.semmod <- function (object, file = "", text, ...) {
	regularizePaths <- function(p) sub(" *-*> *", "->", sub(" *<-* *", "<-", p))
	delete.model.element <- function(delete.text, old.model, type = "path") {
		type <- match.arg(type, c("path", "variable", "coefficient"))
		if (type == "path") delete.text <- regularizePaths(delete.text)
		col.index <- c(path = 1, variable = 1, coefficient = 2)[type]
		delete.text <- strip.white(delete.text)
		old.model <- old.model[-grep(delete.text, strip.white(old.model[, col.index])),]
		class(old.model) <- "semmod"
		return(old.model)
	}
	object[, 1] <- regularizePaths(object[, 1])
	modmat <- if (!missing(text)) scan(text = text, what = list(change = "", var1 = "", 
	                                         var2 = "", var3 = "", var4 = ""), sep = ",", strip.white = TRUE, 
	                                        comment.char = "#", fill = TRUE)
    else scan(file = file, what = list(change = "", var1 = "", 
					var2 = "", var3 = "", var4 = ""), sep = ",", strip.white = TRUE, 
			comment.char = "#", fill = TRUE)
	modmat <- cbind(modmat$change, modmat$var1, modmat$var2, 
			modmat$var3)
	if ("delete" %in% modmat[, 1]) {
		deletemat <- modmat[which(modmat[, 1] == "delete"), 2:3, drop=FALSE]
		deletemat[which(deletemat[, 2] == ""), 2] <- "path"
		for (i in 1:nrow(deletemat)) 
			object <- delete.model.element(deletemat[i, 1], object, deletemat[i, 2])
	}
	if ("add" %in% modmat[, 1]) {
		addmat <- modmat[which(modmat[, 1] == "add"), 2:4, drop=FALSE]
		addmat[addmat[, 3] == "", 3] <- NA
		addmat[addmat[, 2] == "", 2] <- NA
		class(addmat) <- "semmod"
		object <- combineModels(object, addmat, warn=FALSE)
	}
	if ("replace" %in% modmat[, 1]) {
		submat <- modmat[which(modmat[, 1] == "replace"), 2:3, drop=FALSE]
		for (i in 1:nrow(submat)){
			object[, 1:2] <- mapply(function(x, y) gsub(x, y, object[, 1:2]), submat[i, 1], submat[i, 2])
		}
	}
	if ("fix" %in% modmat[, 1]){
		fixmat <- modmat[which(modmat[, 1] == "fix"), 2:3, drop=FALSE]
		fixmat[, 1] <- regularizePaths(fixmat[, 1])
		for (i in 1:nrow(fixmat)){
			object[which(object[, 1] == fixmat[i, 1]), 2:3] <- c(NA, fixmat[i, 2])
		}
	}
	if ("free" %in% modmat[, 1]){
		freemat <- modmat[which(modmat[, 1] == "free"), 2:4, drop=FALSE]
		freemat[, 1] <- regularizePaths(freemat[, 1])
		for (i in 1:nrow(freemat)){
			object[which(object[, 1] == freemat[i, 1]), 2:3] <- freemat[i, 2:3]
		}
	}
	removeRedundantPaths(object)
}

# the next function is adapted from code contributed by Michael Friendly

edit.semmod <- function(name, ...) {
	name <- unclass(name)
	colnames(name) <- c("Path", "Parameter", "StartValue")
	name <- edit(name)
	class(name) <- "semmod"
	name
}
