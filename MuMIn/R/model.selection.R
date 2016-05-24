`coefTable.model.selection` <-
function (model, ...) {
	rval <- attr(model, "coefTables")
	names(rval) <- rownames(model)
	rval
}

`coef.model.selection` <-
function (object, ...) {
	ct <- attr(object, "coefTables")
	n <- length(ct)
	allcf <- unique(unlist(lapply(ct, rownames)))
	rval <- matrix(NA_real_, nrow = n, ncol = length(allcf),
		dimnames = list(rownames(object), allcf))
	for(i in seq_len(n))
		rval[i, match(rownames(ct[[i]]), allcf)] <- ct[[i]][, 1L]
	rval
}

`coeffs.model.selection` <-
function (model) coef.model.selection(model)

`coefArray` <- function(object) {
	coefNames <- fixCoefNames(unique(unlist(lapply(object, rownames),
		use.names = FALSE)))
	nCoef <- length(coefNames)
	nModels <- length(object)
	rval <- array(NA_real_, dim = c(nModels, 3L, nCoef),
		dimnames = list(names(object), c("Estimate", "Std. Error", "df"), coefNames))
	for(i in seq_along(object)) {
		z <- object[[i]]
		rval[i, seq_len(ncol(z)), ] <- t(z[match(coefNames, fixCoefNames(rownames(z))), ])
	}
	rval
}

`getCall.model.selection`  <-
function (x, i = NULL, ...) {
	if(is.null(i)) return(attr(x, "call", exact = TRUE))
	if(length(i) == 1L) return(attr(x, "model.calls", exact = TRUE)[[i]])
	return(attr(x, "model.calls", exact = TRUE)[i])
}

getModelClass <-
function(x) {
	if(inherits(x, "model.selection")) {
		if(!is.null(attr(x, "global"))) return(class(attr(x, "global"))[1L])
		if("class" %in% colnames(x)) return(as.character(x[, "class"]))
		if(!is.null(attr(x, "model.class"))) return(attr(x, "model.class"))
	}
	return(NULL)
}

`print.model.selection` <-
function(x, abbrev.names = TRUE, warnings = getOption("warn") != -1L, ...) {
	origx <- x
	class(x) <- "data.frame"
	xterms <- attr(origx, "terms")
	if(is.null(xterms) || !all(xterms %in% colnames(x)[seq_along(xterms)])) {
		print.data.frame(x, ...)
	} else {
		if(abbrev.names) xterms <- abbreviateTerms(xterms, 6L, 3L, deflate = TRUE)
		colnames(x)[seq_along(xterms)] <- xterms
		globcl <- attr(origx, "global.call")
		if(!is.null(globcl)) {
			cat("Global model call: ")
			print(globcl)
			cat("---\n")
			random.terms <- attr(getAllTerms(attr(origx, "global")), "random.terms")
			if(!is.null(random.terms)) random.terms <- list(random.terms)
		} else random.terms <- attr(origx, "random.terms")
		cat("Model selection table \n")
		
		dig <- c(terms = NA, varying = NA, extra = NA, df = 0L, loglik = 3L, ic = 1L, delta = 2L,
				 weight = 3L)
		column.types <- attr(origx, "column.types")
		#stopifnot(names(dig) == levels(column.types)) ## DEBUG
		
		decprint <- dig[column.types[colnames(x)]]

		i <- vapply(x, is.numeric, FALSE) & is.na(decprint)
		x[, i] <- signif(x[, i], 4L)
		k <- which(!is.na(decprint))
		for(i in k) x[, i] <- round(x[, i], digits = decprint[i])
			
		vLegend <- NULL
		if(abbrev.names) {
			vCols <- type2colname(column.types, "varying")
			vCols <- vCols[(vCols %in% colnames(x)) & !(vCols %in% c("class"))]
			vlen <- nchar(vCols)
			vLegend <- vector(length(vCols), mode = "list")
			names(vLegend) <- vCols
			if(!is.null(vCols)) {
				for(i in vCols) {
					if(!is.factor(x[, i])) next
					lev <- levels(x[, i])
					lev <- lev[!(lev %in% c("", "NULL"))]
					shlev <- abbreviateTerms(lev, nchar(i), deflate = TRUE)
					x[, i] <- factor(x[, i], levels = lev, labels = shlev)
					if(any(j <- shlev != lev)) vLegend[[i]] <-
						paste(shlev[j], "=", sQuote(lev[j]))
				}
				vLegend <- vLegend[!vapply(vLegend, is.null, TRUE)]
			}
		}

		uqran <- unique(unlist(random.terms, use.names = FALSE))
		abbran <- abbreviateTerms(gsub("1 | ", "", uqran, fixed = TRUE), 1L,
			deflate = TRUE)
		colran <- vapply(random.terms, function(s) paste(abbran[match(s, uqran)],
			collapse = "+"), "")

		if(addrandcol <- length(unique(colran)) > 1L) {
			k <- which(colnames(x) == "df")[1L]
			x <- cbind(x[, 1L:(k - 1L)], random = colran, x[, k:ncol(x)])
		}

		if(nrow(x) == 0L) {
			print.default(colnames(x), quote = FALSE)
			cat("<0 rows>", "\n")
		} else
			print.default(as.matrix(x)[, !vapply(x, function(y) all(is.na(y)), FALSE),
			drop = FALSE], na.print = "", quote = FALSE, right = TRUE)

		if(abbrev.names && length(vLegend)) {
			cat("Abbreviations:", sep = "\n")
			for(i in names(vLegend)) {
				cat(vLegend[[i]], sep = ", ", fill = TRUE, labels =
					c(paste0(i, ":"), rep(paste(rep(" ", nchar(i) + 1L),
					collapse = ""), length(vLegend[[i]]) - 1L)))
			}
		}
		
		cat("Models ranked by", asChar(attr(attr(origx, 'rank'), "call")), "\n")
		if(!is.null(random.terms)) {
			if(addrandcol) {
				cat("Random terms: \n")
				cat(paste(abbran, "=", sQuote(uqran)), sep = "\n")
			} else {
				cat("Random terms (all models): \n")
				cat(paste(sQuote(uqran)), sep = ", ")
				cat("\n")
			}
		}
		if (warnings && !is.null(attr(origx, "warnings"))) {
			cat("\n"); print.warnings(attr(origx, "warnings"))
		}
	}
	invisible(origx)
}

`update.model.selection` <- function (object, global.model, ..., evaluate = TRUE) {
    cl <- attr(object, "call")
    if (is.null(cl)) stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...

	if(!missing(global.model))
		extras <- c(list(global.model = substitute(global.model)), extras)

    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(cl)))
        for (a in names(extras)[existing]) cl[a] <- extras[a]
        if (any(!existing)) {
            cl <- c(as.list(cl), extras[!existing])
            cl <- as.call(cl)
        }
    }
    if (evaluate) eval.parent(cl) else cl
}

`logLik.model.selection` <- function (object, ...) {
	nobs <- attr(object, "nobs")
	n <- nrow(object)
	rval <- vector(n, mode = "list")
	for(i in 1L:n) rval[[i]] <-
		structure(object[i, "logLik"], df = object[i, "df"], nobs = nobs,
			class = "logLik")
	rval
}

`family.model.selection` <-
function (object, ...) {
	if(!is.null(attr(object, "global"))) {
		model.calls <- attr(object, "model.calls")
		if(!is.null(model.calls[[1L]][["family"]])) {
			fam <- lapply(model.calls, "[[", "family")
			#fam1 <- unique(fam)
			rval <- lapply(unique(fam), eval)[
				as.integer(as.factor(vapply(fam, asChar, "")))
				]
			names(rval) <- rownames(object)
			## WTF?
			#index <- split(seq_along(fam), vapply(fam, asChar, ""))
			#for(i in seq_along(fam1)) fam1[[i]] <- list(family = eval(fam1[[i]]), index = index[[i]])
			#fam <- family(dd1)
			#index <- lapply(fam, "[[", "index")
			#rval <- rep(lapply(fam, "[[", "family"), vapply(index, length, 1L))[order(unlist(index))]
			return(rval)
		} else return(family(attr(object, "global")))
	} else {
		attr(object, "model.family")
	}
}

`nobs.model.selection` <-
function (object, ...)
attr(object, "nobs")

## internal: translate column type to column indices
type2col <-
function (x, type) {
    if (inherits(x, "model.selection")) 
        x <- attr(x, "column.types")
	k <- match(x, type, nomatch = 0L)
	i <- k != 0
	which(i)[order(k[i])]
}

## internal: translate column type to column names
type2colname <-
function(x, type)
names(x)[type2col(x, type)] 


`item<-` <- function(x, name, i, value)
`[<-.data.frame`(x, i, name, value)

`item` <- function(x, name, i, ...)
`[.data.frame`(x, i, name, ...)

`itemByType` <- function(x, type, i, ...) 
`[.data.frame`(x, i, type2col(x, type), ...)

`itemByType<-` <- function(x, type, i, value)
`[<-.data.frame`(x, i, type2col(x, type), value)
