`row.names<-.model.selection` <-
function (x, value)  {
	oldnames <- dimnames(x)[[1L]]
	x <- NextMethod()
	newnames <- dimnames(x)[[1L]]
	rowattrib <- c("model.calls", "coefTables", "random.terms", "order",
		if(!is.null(attr(x, "modelList"))) "modelList")
	for(i in rowattrib) if(!is.null(attr(x, i))) names(attr(x, i)) <- newnames
	x
}

`names<-.model.selection` <-
function (x, value) {
	oldnames <- names(x)
	if(any(attr(x, "column.types")[oldnames[oldnames != value]] %in%
	   c('df', 'loglik', 'ic', 'delta', 'weight', 'terms'))) {
		class(x) <- "data.frame"
		attributes(x)[-match(names(attributes(x)), c("names", "row.names", "class"),
							 nomatch = 0)] <- NULL
	}
	NextMethod()
}

subset_model_selection <-
function(x, attrib, modif = NULL, rowchange = TRUE) {
	excludeattr <- c("names", "row.names", "class")
	column.types <- attrib[["column.types"]]
	keepattr <- names(attrib)[!(names(attrib) %in% excludeattr)]
	.setattr <- function(x, newattr = NULL, which = keepattr) {
		attributes(x)[which] <- if(is.null(newattr)) NULL else newattr[which]
		x
	}

	if(inherits(x, "model.selection")) {
		protectedcoltypes <- c("df", "loglik", "ic", "delta", "weight", "terms")

		if(!is.null(modif) && modif %in% type2colname(column.types, protectedcoltypes)) {
			class(x) <- "data.frame"
			return(.setattr(x))
		} else {
			s <- dimnames(x)[[2L]]
			k <- match(names(column.types), colnames(x), nomatch = 0L)
			if(any(column.types[k == 0L] %in% protectedcoltypes)) {
				class(x) <- "data.frame"
				return(.setattr(x))
			} else {
				if(any(column.types[k == 0L] %in% c("varying", "extra"))) {
					column.types <- column.types[k != 0L]
					attrib[["column.types"]] <- column.types
				}
			}
		}
		oldrownames <- attrib[['row.names']]
		newrownames <- dimnames(x)[[1L]]
		if(rowchange && (length(oldrownames) != length(newrownames) ||
						 any(oldrownames != newrownames))) {
			rowattrib <- c("model.calls", "coefTables", "random.terms", "order",
			   if(!is.null(attr(x, "modelList")))"modelList")
			k <- match(newrownames, oldrownames)
			attrib[rowattrib] <- lapply(attrib[rowattrib], `[`, k)
		}
		x <- .setattr(x, attrib)
		if(!is.null(warningList <- attrib$warnings))
			attr(x, "warnings") <- warningList[sapply(warningList, attr, "id")
											   %in% newrownames]
	} else {
		return(.setattr(x))
	}
	x
}


`[<-.model.selection` <-
function (x, i, j, value)  {
	if (missing(j)) j <- TRUE
	subset_model_selection(NextMethod("[<-"),
		attributes(x), if(is.character(j)) j else colnames(x)[j])
}


`[[<-.model.selection` <-
function (x, i, j, value)  {
	subset_model_selection(NextMethod(),
		attributes(x), {
			if(missing(j)) j <- i
			if(is.character(j)) j else colnames(x)[j]
		}, rowchange = FALSE)
}

`$<-.model.selection` <-
function (x, name, value) {
	subset_model_selection(NextMethod("$<-"), attributes(x), name,
		rowchange = FALSE)
}

`[.model.selection` <-
function (x, i, j, recalc.weights = TRUE, recalc.delta = FALSE, ...) {
	x <- subset_model_selection(item(x, j, i, ...),	origattrib <- attributes(x))
	if(inherits(x, "model.selection")) {
		ic <- itemByType(x, "ic")
		if(recalc.weights) itemByType(x, "weight") <- Weights(ic)
		if(recalc.delta) itemByType(x, "delta") <- ic - min(ic)
	} else {
		k <- type2colname(origattrib$column.types, c("weight", "delta"))
		hasdeltaweight <- k %in% colnames(x)
		recalc <- c(if(recalc.delta && hasdeltaweight[2L]) "delta",
					if(recalc.weights && hasdeltaweight[1L]) "weights")
		if(!is.null(recalc)) cry(, "cannot recalculate %s on an incomplete object",
					prettyEnumStr(recalc), warn = TRUE)
	}
	x
}


`[[.model.selection` <-
function (x, ..., exact = TRUE) {
	`[[.data.frame`(x, ..., exact = exact)
}


subset_rework <-
function(subset, object, objectname = substitute(object)) {

	if(!is.language(subset) && is.na(subset)) return(TRUE)

	subset <- exprapply0(exprapply0(exprapply0(subset, "dc", .sub_dc_has),
		c("{", "Term"), .sub_Term),
		    "has", .sub_has)
			
	objectname <- as.name(objectname)

	subset <- exprApply(subset, names(object), symbols = TRUE,
		function(x, v, cl, parent) {
			if(is.call(parent) && any(parent[[1L]] == c("I", "$", "@")))
				return(x)
			if(length(x) == 1L) {
				cl[[3L]] <- match(asChar(x), v)
				return(cl)
			}
			x
		}, v = names(object), call("[[", objectname, 0L))
	
	subset <- exprApply(subset, "I", function(x) x[[2L]])
	subset <- subst(subset, . = objectname)
	subset
}

subset_eval <-
function(subset, x, envir) {
	eval(subset_rework(subset, x, "tmp_data_"),
		list(tmp_data_ = x, dc = .subset_vdc),
		envir
		)
}


`subset.model.selection` <-
function(x, subset, select, recalc.weights = TRUE, recalc.delta = FALSE, ...) {
	if(missing(subset) && missing(select)) return(x)
	
	#ss <- eval(subset_rework(substitute(subset), x, "tmpdat"),
	#	list(tmpdat = x, .subset_vdc = .subset_vdc),
	#	parent.frame()
	#	)
	
	# eval.parent(subset_rework(substitute(subset), x, substitute(x)))
	
	return(`[.model.selection`(x,
		subset_eval(substitute(subset), x, parent.frame()),
		recalc.weights = recalc.weights,
		recalc.delta = recalc.delta, ...))
}
