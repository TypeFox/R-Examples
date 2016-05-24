xtable.summary.averaging <-
function (x, caption = NULL, label = NULL, align = NULL, digits = NULL, 
    display = NULL, coefType = c("full", "subset"), ...) {
	coefType <- match.arg(coefType)
    x <- as.data.frame(x[[switch(coefType, full = "coefmat.full",
		subset = "coefmat.subset")]])
	has.ase <- all(!is.na(x[, 3L]))
	if(!has.ase) x <- x[, -3L]
	getFrom("xtable", "xtable")(x, caption = caption, label = label,
		   align = if(is.null(align)) rep("r", ncol(x) + 1L) else align,
		   digits = if(is.null(digits)) c(0, 4, 4, if(has.ase) 4, 2, 4) else digits, 
		display = if(is.null(display)) c("s", "f", "f", if(has.ase) "f", "f", "f") else display,
		...)
    return(x)
}

xtable.averaging <- 
function (x, caption = NULL, label = NULL, align = NULL, digits = NULL, 
    display = NULL, coefType, ...) {
    return(xtable.summary.averaging(summary(x, ...), caption = caption, label = label, 
        align = align, digits = digits, display = display, coefType = coefType))
}

xtable.model.selection <-
function (x, caption = NULL, label = NULL, align = NULL, digits = NULL,
		  display = NULL, ...) {
	column.types <- attr(x, "column.types")
	x <- as.data.frame(x)
	vclass <- vapply(x, function(v) {
		if(is.integer(v)) return("integer")
		if(is.factor(v) || is.character(v)) return("character")
		if(is.numeric(v)) return("real")
		"other"
	}, "")
	if(is.null(align)) align <- c("r", ifelse(vclass == "character", "c", "r"))
	
	dig <- c(terms = NA, varying = NA, extra = NA, df = 0L, loglik = 1L,
				 ic = 1L, delta = 1L, weight = 2L)
	decprint <- dig[column.types[colnames(x)]]
	decprint[is.na(decprint)] <- 2L
	display <- character(ncol(x))
	display[vclass == "character"] <- "s"
	display[vclass == "real"] <- "f"
	display[vclass == "integer"] <- "d"
	display[vclass == "other"] <- "s"
	getFrom("xtable", "xtable")(x, caption = caption, label = label, align = align, digits = c(NA, decprint), 
		display = c("s", display), ...)
}





