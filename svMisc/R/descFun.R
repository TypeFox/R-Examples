descFun <- function (fun, package, lib.loc = NULL)
{
	if (!length(fun)) return("")
	fun <- as.character(fun)
	l <- length(fun)
	if (missing(package) || is.null(package)) package <- ""
	package <- rep(package, length.out = l)
	
	## Create a vector of results
	res <- rep("", l)
	
	## Collect help for each function
	for (i in 1:l) {
		## Get location of the help file
		## We cannot just call help normally because otherwise it thinks
		## we are looking for package "package" so we create a call and eval it
		help.call <- call("help", fun[i], lib.loc = lib.loc, help_type = "text")
		if (package[i] != "") help.call[["package"]] <- package[i]
		file <- eval(help.call)
		file <- as.character(file)
		if (length(file) > 0) {
			## Read the Rd file and get the title section out of it
			Rdoc <- getNamespace("utils")$.getHelpFile(file[1L])
			## Look for the \title tag
			j <- 0
			for (j in seq_along(Rdoc))
				if (attr(Rdoc[[j]], "Rd_tag") == "\\title") break
			if (j > 0) {
				desc <- as.character(Rdoc[[j]][[1]])
				desc <- sub("^[ \t]+", "", desc)
				desc <- sub("[ \t]+$", "", desc)
				res[i] <- desc
			}
		}
	}
	res
}

descArgs <- function (fun, args = NULL, package = NULL, lib.loc = NULL)
{	
	## We cannot just call help normally because otherwise it thinks
	## we are looking for package "package" so we create a call and eval it
	help.call <- call("help", fun, lib.loc = lib.loc, help_type = "text")
	if (!is.null(package)) help.call[["package"]] <- package
	file <- eval(help.call)
	
	## This is borrowed from utils::print.help_files_with_topic
	path <- dirname(file)
    dirpath <- dirname(path)
    pkgname <- basename(dirpath)
    RdDB <- file.path(path, pkgname)
    
    if (!file.exists(paste(RdDB, "rdx", sep=".")))
    	return(character(length(args)))
    
	rd <- getNamespace("tools")$fetchRdDB(RdDB, basename(file))
    
    ## This is not exported from tools
    RdTags <- function (Rd) {
    	res <- sapply(Rd, attr, "Rd_tag")
    	if (!length(res)) res <- character(0)
    	return(res)
    }
    tags <- gsub("\\", "", RdTags(rd), fixed = TRUE) 
    
    if (!any(tags == "arguments")) return(character(length(args)))
    
    arguments <- rd[[which(tags == "arguments")[1]]]
    items <- arguments[RdTags(arguments) == "\\item"]
    descriptions <- do.call(rbind, lapply(items, function (item) {
		names <- try(strsplit(item[[1]][[1]], "\\s*,\\s*", perl = TRUE)[[1]],
			silent = TRUE)
		if (inherits(names, "try-error")) {
			## This happens with the "..." argument
			names <- "..."
		}
    	content <- paste(rapply(item[-1], as.character), collapse = "")
    	cbind(names, rep.int(content, length(names)))
    }))
    
    if (is.null(args)) {
    	structure(descriptions[, 2], names = descriptions[, 1])
    } else {
    	sapply(args, function (a) {
    		if (a %in% descriptions[, 1]) {
    			descriptions[which(descriptions[, 1] == a)[1] , 2]
    		} else ""
    	})
    }
}

.descData <- function (data, columns, package = NULL, lib.loc = NULL)
	character(length(columns))

.descSlots <- function (object, slots, package = NULL, lib.loc = NULL)
	character(length(slots))

.descSquare <- function (completions, package = NULL)
	character(length(completions))
