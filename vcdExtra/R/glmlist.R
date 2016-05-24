# glmlist - make a glmlist object containing a list of fitted glm objects with their names

# borrowing code from Hmisc::llist
glmlist <- function(...) {
    args  <- list(...);
    lname <- names(args)
    name <- vname <- as.character(sys.call())[-1]
    for (i in 1:length(args)) {
        vname[i] <- if (length(lname) && lname[i] != "") 
            lname[i]
        else name[i]
    }
    names(args) <- vname[1:length(args)]
    is.glm <- unlist(lapply(args, function(x) inherits(x, "glm")))
    if (!all(is.glm)) {
    	warning("Objects ", paste(vname[!is.glm], collapse=', '), 
              " removed because they are not glm objects")
    	args <- args[is.glm]
    }
    class(args) <- "glmlist"
    return(args);
}

# loglmlist - do the same for loglm objects

loglmlist <- function(...) {
    args  <- list(...);
    lname <- names(args)
    name <- vname <- as.character(sys.call())[-1]
    for (i in 1:length(args)) {
        vname[i] <- if (length(lname) && lname[i] != "") 
            lname[i]
        else name[i]
    }
    names(args) <- vname[1:length(args)]
    is.loglm <- unlist(lapply(args, function(x) inherits(x, "loglm")))
    if (!all(is.loglm)) {
    	warning("Objects ", paste(vname[!is.loglm], collapse=', '), 
              " removed because they are not loglm objects")
    	args <- args[is.loglm]
    }
    class(args) <- "loglmlist"
    return(args);
}

# generic version: named list

nlist <- function(...) {
    args  <- list(...);
    lname <- names(args)
    name <- vname <- as.character(sys.call())[-1]
    for (i in 1:length(args)) {
        vname[i] <- if (length(lname) && lname[i] != "") 
            lname[i]
        else name[i]
    }
    names(args) <- vname[1:length(args)]
    return(args);
}

# coeficient method for a glmlist (from John Fox, r-help, 10-28-2014)

coef.glmlist <- function(object, result=c("list", "matrix", "data.frame"),
		...){
	result <- match.arg(result)
	coefs <- lapply(object, coef)
	if (result == "list") return(coefs)
	coef.names <- unique(unlist(lapply(coefs, names)))
	n.mods <- length(object)
	coef.matrix <- matrix(NA, length(coef.names), n.mods)
	rownames(coef.matrix) <- coef.names
	colnames(coef.matrix) <- names(object)
	for (i in 1:n.mods){
		coef <- coef(object[[i]])
		coef.matrix[names(coef), i] <- coef
	}
	if (result == "matrix") return(coef.matrix)
	as.data.frame(coef.matrix)
}

