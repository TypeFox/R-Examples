objList <- function (id = "default", envir = .GlobalEnv, object = NULL,
all.names = FALSE, pattern = "", group = "", all.info = FALSE, sep = "\t",
path = NULL, compare = TRUE, ...)
{
	## Make sure that id is character
	id <- as.character(id)[1]
	if (id == "") id <- "default"
	ename <- NA

	## Format envir as character (use only first item provided!)
	if (!is.environment(envir)){
		if (is.numeric(envir) && envir > 0)
			envir <- search()[envir]

		if (is.character(envir)) {
			ename <- envir
			envir <- tryCatch(as.environment(envir), error = function(e) NULL)
			if (is.null(envir) || inherits(envir, "error")) {
				envir <- NULL
				ename <- ""
			}
		}
	}

	# base and .GlobalEnv do not have name attribute
	if (!is.null(attr(envir, "name"))) ename <- attr(envir, "name")
	else if (is.na(ename)) ename <- deparse(substitute(envir))
	if (ename %in% c("baseenv()", ".BaseNamespaceEnv"))
		ename <- "package:base"


	## Object to return in case of empty data
	#Nothing <- data.frame(Envir = character(0), Name = character(0),
	#	Dims = character(0), Group = character(0), Class = character(0),
	#	Recursive = logical(0), stringsAsFactors = FALSE)
	#if (!isTRUE(all.info)) Nothing <- Nothing[, -1]
	#attr(Nothing, "all.info") <- all.info
	#attr(Nothing, "envir") <- ename
	#attr(Nothing, "object") <- object
	#class(Nothing) <- c("objList", "data.frame")

	# This is ~15x faster:
	Nothing <- structure(list(Name = character(0),
		Dims = character(0), Group = character(0), Class = character(0),
		Recursive = logical(0), stringsAsFactors = FALSE),
			class=c("objList", "data.frame"),
			all.info= all.info, envir=ename, object=object
		)
		if (isTRUE(all.info)) Nothing <- cbind(Envir = character(0), Nothing)


	if (is.null(envir)) return(Nothing)

	if (!missing(object) && is.character(object) && object != "") {
		res <- .lsObj(envir = envir, objname = object)
	} else {
		## Get the list of objects in this environment
		Items <- ls(envir = envir, all.names = all.names, pattern = pattern)
		if (length(Items) > 0) {
			## Get characteristics of all objects
			describe <- function (name, all.info = FALSE) {
				## Get a vector with five items:
				## Name, Dims, Group, Class and Recursive
				obj <- envir[[name]]
				res <- c(
					Name = name,
					Dims = if (is.null(Dim <- dim(obj))) length(obj) else
						paste(Dim, collapse = "x"),
					Group = mode(obj),
					Class = class(obj)[1],
					Recursive = is.recursive(obj) || mode(obj) == "S4"
				)
				return(res)
			}
			res <- data.frame(t(sapply(Items, describe, all.info = all.info)),
				stringsAsFactors = FALSE)
	
			# Quote non-syntactic names
			nsx <- res$Name != make.names(res$Name)
			res$Full.name[!nsx] <- res$Name[!nsx]
			res$Full.name[nsx] <- paste("`", res$Name[nsx], "`", sep = "")
			res <- res[, c(1, 6, 2:5)]
		} else res <- Nothing
	
		## No, because if rm(list = ls()), we must reactualize the objects
		## browser anyway
		#if (NROW(res) == 0) return(Nothing)
	
		if (isTRUE(all.info)) res <- cbind(Envir = ename, res)
	
		vMode <- Groups <- res$Group
		vClass <- res$Class
	
		## Recalculate groups into meaningful ones for the object explorer
		## 1) Correspondance of typeof() and group depicted in the browser
		Groups[Groups %in% c("name", "environment", "promise", "language", "char",
			"...", "any", "(", "call", "expression", "bytecode", "weakref",
			"externalptr")] <- "language"
	
		Groups[Groups == "pairlist"] <- "list"
	
		## 2) All Groups not being language, function or S4 whose class is
		##    different than typeof are flagged as S3 objects
		Groups[!(Groups %in% c("language", "function", "S4")) &
			vMode != vClass] <- "S3"
	
		## 3) Integers of class factor become factor in group
		Groups[vClass == "factor"] <- "factor"
	
		## 4) Objects of class 'data.frame' are also group 'data.frame'
		Groups[vClass == "data.frame"] <- "data.frame"
	
		## 5) Objects of class 'Date' or 'POSIXt' are of group 'DateTime'
		Groups[vClass == "Date" | vClass == "POSIXt"] <- "DateTime"
	
		## Reaffect groups
		res$Group <- Groups
	
		## Possibly filter according to group
		if (!is.null(group) && group != "")
			res <- res[Groups == group, ]
	}

	## Determine if it is required to refresh something
	Changed <- TRUE
	if (isTRUE(compare)) {
		allList <- getTemp(".guiObjListCache", default = list())

		if (identical(res, allList[[id]])) Changed <- FALSE else {
			## Keep a copy of the last version in SciViews:TempEnv
			allList[[id]] <- res
			assignTemp(".guiObjListCache", allList)
		}
	}

	## Create the 'objList' object
	attr(res, "all.info") <- all.info
	attr(res, "envir") <- ename
	attr(res, "object") <- object
	attr(res, "changed") <- Changed
	attr(res, "class") <- c("objList", "data.frame")

	if (is.null(path)) {  # Return results or "" if not changed
		return(if (Changed) res else Nothing)
	} else if (Changed) {  # Write to files in this path
		return(write.objList(res, path = path, sep = sep, ...))
	} else {
		return(Nothing)  # Not changed
	}
}

write.objList <- function (x, path, sep = "\t", ...)
{
	id <- attr(x, "id")
	ListF <- file.path(path, sprintf("List_%s.txt", id))
	ParsF <- file.path(path, sprintf("Pars_%s.txt", id))

	write.table(as.data.frame(x), row.names = FALSE, col.names = FALSE,
		sep = sep, quote = FALSE, file = ListF)

	## Write also in the Pars_<id>.txt file in the same directory
	cat(sprintf("envir=%s\nall.names=%s\npattern=%s\ngroup=%s",
		attr(x, "envir"), attr(x, "all.names"), attr(x, "pattern"),
		attr(x, "group")), file = ParsF, append = FALSE)

	invisible(ListF)
}

print.objList <- function (x, sep = NA, eol = "\n",
header = !attr(x, "all.info"), raw.output = !is.na(sep), ...)
{
	if (!inherits(x, "objList"))
		stop("x must be an 'objList' object")

	empty <- NROW(x) == 0

	if (!raw.output)
		cat(if (empty) "An empty objects list\n" else "Objects list:\n")

	if (header) {
		header.fmt <- if (raw.output) "Env=%s\nObj=%s\n" else
			"\tEnvironment: %s\n\tObject: %s\n"

		objname <- if (is.null(attr(x, "object"))) {
			if (raw.output) "" else "<All>"
		} else attr(x, "object")

		cat(sprintf(header.fmt,  attr(x, "envir"), objname))
	}

	if (!empty) {
		if (is.na(sep)) {
			print(as.data.frame(x))
		} else if (!is.null(nrow(x)) && nrow(x) > 0) {
			write.table(x, row.names = FALSE, col.names = FALSE, sep = sep,
				eol = eol, quote = FALSE)
		}
	}
	invisible(x)
}

## Called by objList() when object is provided
.lsObj <- function (objname, envir, ...)
{
	obj <- try(eval(parse(text = objname), envir = as.environment(envir)),
		silent = TRUE)
	if (inherits(obj, "try-error")) return(NULL)

	if (is.environment(obj)) obj <- as.list(obj)

	if (mode(obj) == "S4") {
		ret <- .lsObjS4(obj, objname)
	} else if (is.function(obj)) {
		ret <- .lsObjFunction(obj, objname)
	} else {  # S3
		if (!(mode(obj) %in% c("list", "pairlist")) || length(obj) == 0)
			return(NULL)

		itemnames <- fullnames <- names(obj)
		if (is.null(itemnames)) {
			itemnames <- seq_along(obj)
			fullnames <- paste(objname, "[[", itemnames, "]]", sep = "")
		} else {
			w.names <- itemnames != ""
			.names <- itemnames[w.names]
			nsx <- .names != make.names(.names)  # Non-syntactic names
			.names[nsx] <- paste("`", .names[nsx], "`", sep = "")
			fullnames[w.names] <- paste (objname, "$", .names, sep = "")
			fullnames[!w.names] <- paste(objname, "[[",
				seq_along(itemnames)[!w.names], "]]", sep = "")
		}

		ret <- t(sapply(seq_along(obj), function (i) .objDescr(obj[[i]])))

		ret <- data.frame(itemnames, fullnames, ret, stringsAsFactors = FALSE)
	}
	if (!is.null(ret))
		names(ret) <- c("Name", "Full.name", "Dims/default", "Group", "Class",
			"Recursive")
	ret
}

# Called by .lsObj for functions
.lsObjFunction <- function (obj, objname = deparse(substitute(obj)))
{
	## formals(obj) returns NULL if only arg is ..., try: formals(expression)
	obj <- formals(args(obj))
	objname <- paste("formals(args(", objname, "))", sep = "")

	if(length(obj) == 0) return(NULL)

	itemnames <- fullnames <- names(obj)
	nsx <- itemnames != make.names(itemnames) # non-syntactic names
	itemnames[nsx] <- paste("`", itemnames[nsx], "`", sep = "")
	fullnames <- paste(objname, "$", itemnames, sep = "")

	ret <- t(sapply (seq_along(obj), function (i) {
		x <- obj[[i]]
		lang <- is.language(obj[[i]])
		o.class <- class(obj[[i]])[1]
		o.mode <- mode(obj[[i]])

		d <- deparse(obj[[i]])
		if (lang && o.class == "name") {
			o.class <- ""
			o.mode <- ""
		}

		ret <- c(paste(d, collapse = "x"), o.class,	o.mode, FALSE)
		return(ret)
	}))

	ret <- data.frame(itemnames, fullnames, ret, stringsAsFactors = FALSE)
	ret
}

## Called by .lsObj in S4 case
.lsObjS4 <- function (obj, objname = deparse(substitute(obj)))
{
	itemnames <- fullnames <- slotNames(obj)
	nsx <- itemnames != make.names(itemnames)
	itemnames[nsx] <- paste("`", itemnames[nsx], "`", sep = "")
	fullnames <- paste(objname, "@", itemnames, sep = "")

	ret <- t(sapply(itemnames, function (i) .objDescr(slot(obj, i))))

	ret <- data.frame(itemnames, fullnames, ret, stringsAsFactors = FALSE)
	ret
}

## Returns a *character* vector with elements: dims, mode, class, rec(ursive)
.objDescr <- function (x) {
	d <- dim(x)
	if (is.null(d)) d <- length(x)

	c(dims = paste(d, collapse = "x"),
		mode = mode(x), class = class(x)[1],
		rec = mode(x) == "S4" || is.function(x) ||
		(is.recursive(x) && !is.language(x) && sum(d) != 0))
}
