## TODO: complex => character + how to restore complex numbers with attributes = TRUE?
## TODO: check dates, and manage other dates than Date!
## TODO: convert functions, expressions into string, and how to include JS code? or R code?
## TODO: allow for special characters \b, \n, \r, \f, \t, \" in names!
## TODO: environment and proto
toRjson <- function (x, attributes = FALSE) 
{
	## This is derived from dput()
	file <- file()
	on.exit(close(file))
	if (isTRUE(attributes)) {
		opts <- c("showAttributes", "S_compatible")
	} else {
		opts <- "S_compatible"
	}
	
	## Non-named list items are not allowed => make sure we give names to these
	## Also if attributes == FALSE, we use the string representation of factors
	rework <- function (x, attributes = FALSE) {
		if (is.list(x) && length(x)) {
			## Make sure all items have names, and use [[x]] for unnamed items
			i <- paste("[[", 1:length(x), "]]", sep = "")
			n <- names(x)
			if (is.null(n)) {
				n <- i
			} else {
				nonames <- n == ""
				n[nonames] <- i[nonames]
			}
			## Flag names with leading and trailing sequence (unlikely elsewhere)
			n <- paste("@&#&&", n, "&&#&@", sep = "")
			## Change names of x
			names(x) <- n
			## If we don't use attributes, convert factors and Dates to characters
			if (!isTRUE(attributes))
				x <- rapply(x, as.character, classes = c("factor", "Date"),
					how = "replace")
			## Do this recursively
			for (item in names(x))
				x[[item]] <- rework(x[[item]], attributes)
		} else if (!isTRUE(attributes) && inherits(x, c("factor", "Date")))
			x <- as.character(x) 
		## Process also all attributes
		if (isTRUE(attributes)) {
			a <- attributes(x)
			if (!is.null(a)) {
				n <- names(x)
				a$.Names <- NULL
				a$names <- NULL
				na <- names(a)
				if (length(na)) {
					for (item in na)
						a[[item]] <- rework(a[[item]], attributes)
					## Tag attributes names and translate a few special ones
					specials <- c(".Dim", ".Dimnames", ".Tsp", ".Label")
					replace <- c("dim", "dimnames", "tsp", "levels")
					m <- match(na, specials)
					ok <- (!is.na(m) & m)
					na[ok] <- replace[m[ok]]
					names(a) <- paste("@&#&&", na, "&&#&@", sep = "")
				}
				attributes(x) <- a
				names(x) <- n
			}
		}
		return(x)
	}

    ## Is this an S4 object => process each slot separately
	if (isS4(x)) {
		cat('list("Class_" := "', class(x), '"\n', file = file, sep = "")
		for (n in slotNames(x)) {
			cat('    , "', n, '" := ', file = file)
			dput(rework(slot(x, n), attributes), file = file, control = opts)
		}
		cat(")\n", file = file)
		invisible()
	}
	else {
		## Was: .Internal(dput(rework(x, attributes), file, opts))
		dput(rework(x, attributes), file = file, control = opts)
	}
	
	## Now read content from the file
	res <- readLines(file)
	
	## dput() indicates sequences of integers with x:y that JavaScript cannot
	## process... replace these by the equivalent code seq(x, y)
	res <- gsub("(-?[0-9]+):(-?[0-9]+)", "seq(\\1, \\2)", res)
	
	## Convert '.Names = ' into '"names" := '
	res <- gsub(".Names = ", '"names" := ', res, fixed = TRUE)
	## We need to replace special characters
	## TODO: do so only inside `@&#&&...&&#&@`
## TODO: all this does not work!!!
#	res <- gsub('(`@&#&&.*)\b(.*&&#&@`)', '\\1\\\\b\\2', res)
#	res <- gsub('(`@&#&&.*)\t(.*&&#&@`)', '\\1\\\\t\\2', res)
#	res <- gsub('(`@&#&&.*)\n(.*&&#&@`)', '\\1\\\\n\\2', res)
#	res <- gsub('(`@&#&&.*)\f(.*&&#&@`)', '\\1\\\\f\\2', res)
#	res <- gsub('(`@&#&&.*)\r(.*&&#&@`)', '\\1\\\\r\\2', res)
#	res <- gsub('(`@&#&&.*)\"(.*&&#&@`)', '\\1\\\\"\\2', res)
	#res <- gsub('\t', '\\t', res, fixed = TRUE)
	#res <- gsub('\n', '\\n', res, fixed = TRUE)
	#res <- gsub('\f', '\\f', res, fixed = TRUE)
	#res <- gsub('\r', '\\r', res, fixed = TRUE)
	#res <- gsub('\"', '\\"', res, fixed = TRUE)
	## Convert `@&#&& into ", and &&#&@` = into " :=
	res <- gsub('"?`@&#&&', '"', res)
	res <- gsub('&&#&@`\"? =', '" :=', res)
	## Convert "@&#&&[[d]]&&#&@" to "" (non-named items)
	res <- gsub('"@&#&&\\[\\[[1-9][0-9]*]]&&#&@"', '""', res)
	## Convert "@&#&& into " and &&#&@" into "
	res <- gsub('"@&#&&', '"', res, fixed = TRUE)
	res <- gsub('&&#&@"', '"', res, fixed = TRUE)
	## No unnamed items, so, convert 'structure(' into 'list("Data_" := ' 
	res <- gsub("([^a-zA-Z0-9._])structure\\(", '\\1list("Data_" := ', res)
	res <- sub("^structure\\(", 'list("Data_" := ', res)
	## Old code!
	## Convert 'list(' into 'hash('
	#res <- gsub("([^a-zA-Z0-9._])list\\(", "\\1hash(", res)
	#res <- sub("^list\\(", "hash(", res)
	
	## Return  the no quoted results
	noquote(res)
}

evalRjson <- function (rjson) {
	## Our custom list() manages to create list() but also new() or structure() items
	list <- function (Class_, Data_, ...) {
		## If there is a "Class_" argument, create new S4 object
		## Note that "Data_" is ignored in this case!
		if (!missing(Class_)) return(new(Class_, ...))
		## If there is a "_Data_" argument, create a structure
		if (!missing(Data_)) return(structure(Data_, ...))
		## otherwise, create a list
		return(base::list(...))
	}
	
	## To convert RJSON data into a R object, simply evaluate it
	## Note: RJSONp objects will be evaluated correctly too
	## providing the <callback>() exists and can manage a single
	## argument (being the RJSOn object converted to R)
	
	## We need first to convert all ':=' into '='
	eval(parse(text = gsub(":=", "=", rjson, fixed = TRUE)))
}

# Simple JSON for lists containing character strings
listToJson <- function (x) {
	if (!is.list(x) && length(x) == 1L)
		return(encodeString(x, quote = '"'))
	x <- lapply(x, listToJson)
	x <- if (is.list(x) || length(x) > 1L) {
		nms <- names(x)
		if (is.null(nms)) {
			paste('[', paste(x, collapse = ','), ']', sep = "")
		} else {
			paste("{", paste(paste(encodeString(make.unique(nms, sep = '#'),
				quote = '"'), ":", x, sep = ""), collapse = ","), "}", sep = "")
		}
	}
	x
}
