addreg.allref <- function(object, data = environment(object), mono, family, start = NULL) {
    t <- if(missing(data))
        terms(object)
    else terms(object, data = data)
    if(is.null(attr(data, "terms")))
        data <- model.frame(object, data)
    else {
        reorder = match(sapply(attr(t, "variables"), deparse,
            width.cutoff = 500)[-1L], names(data))
        if (any(is.na(reorder)))
            stop("model frame and formula mismatch in addreg.allref()")
        if(!identical(reorder, seq_len(ncol(data))))
            data <- data[, reorder, drop = FALSE]
    }
    int <- attr(t, "response")
	os <- attr(t, "offset")
    
    namD <- names(data)
    for (i in namD) if (is.character(data[[i]]))
        data[[i]] <- factor(data[[i]])
    isF <- vapply(data, function(x) is.factor(x) || is.logical(x), NA)
    isF[int] <- FALSE
	if(!is.null(os)) isF <- isF[-os]
    
	termlist <- attr(t, "term.labels")
    nvar <- length(termlist)
	
	npar <- sum(as.numeric(!isF))
	if (any(isF)) npar <- npar + sum(sapply(data[isF], function(x) nlevels(factor(x)) - 1))
	if (substr(family$family,1,7) == "negbin1") npar <- npar + 1
    
    if (missing(mono)) mono <- rep(FALSE, nvar)
    if (is.null(mono)) mono <- rep(FALSE, nvar)
    monotonic <- rep(FALSE, nvar)
    names(monotonic) <- termlist
    monotonic[mono] <- TRUE
	names(monotonic) <- termlist
    
    allref <- list()
	
	start.new.scale <- NULL
	
	if (!is.null(start)) {
		if (length(start) != npar)
			stop(gettextf("number of values in 'start' is %d should equal %d (number of parameters)",
                length(start), npar), domain = NA)
		start.orig <- start
		start.new.int <- start.orig[1]
		start.new.other <- start.orig[-1]
		if (substr(family$family,1,7) == "negbin1") {
			start.new.scale <- start.new.other[length(start.new.other)]
			start.new.other <- start.new.other[-length(start.new.other)]
		}
		this.start <- 2
	} else {
		start.new.int <- NULL
		start.new.other <- NULL
	}
    
    if (nvar == 0) return(list(allref = allref, terms = t, data = data, monotonic = monotonic, start.new = start))
    for (term in termlist) {
		allref[[term]] <- list()
		term2 <- gsub("`","",term)
        if (!isF[term2]) {
			cont.min <- min(data[[term2]])
			cont.max <- max(data[[term2]])
			if (!is.null(start)) {
				if (family$family != "binomial") {
					if (!monotonic[term]) {
						if (start.orig[this.start] < 0) {
							allref[[term]] <- as.list(2:1)
							start.new.int <- start.new.int + start.orig[this.start] * cont.max
							start.new.other[this.start-1] <- -start.orig[this.start]
						} else {
							allref[[term]] <- as.list(1:2)
							start.new.int <- start.new.int + start.orig[this.start] * cont.min
							start.new.other[this.start-1] <- start.orig[this.start]
						}
					} else {
						allref[[term]][[1]] <- 1
						start.new.int <- start.new.int + start.orig[this.start] * cont.min
						start.new.other[this.start-1] <- start.orig[this.start]
					}
				} else {
                    allref[[term]][[1]] <- 1
                    start.new.int <- start.new.int + start.orig[this.start] * cont.min
                    start.new.other[this.start-1] <- start.orig[this.start]
                }
				this.start <- this.start + 1
			} else {
				allref[[term]][[1]] <- 1
				if(family$family != "binomial" & !monotonic[term]) allref[[term]][[2]] <- 2
			}
			attr(allref[[term]],"type") <- 1
		} else {
			lvls <- levels(factor(data[[term]]))
			nlvls <- nlevels(factor(data[[term]]))
			if (!is.null(start)) {
				start.this <- start.orig[this.start:(this.start + nlvls - 2)]
				if (!monotonic[term]) {
					if (family$family != "binomial") {
						allref[[term]] <- as.list(lvls[order(c(0, start.this))])
						start.new.int <- start.new.int + min(c(0, start.this))
						start.new.other[(this.start - 1):(this.start + nlvls - 3)] <- (c(0, start.this) - min(c(0, start.this)))[-which.min(c(0, start.this))]
						attr(allref[[term]], "type") <- 2
					} else {
						allcombins <- combinat::permn(lvls)
						find.order <- which(sapply(allcombins, function(x) all(x == lvls[order(c(0, start.this))])))
						allref[[term]] <- append(allcombins[find.order], allcombins[-find.order])
						start.new.int <- start.new.int + min(c(0, start.this))
						start.new.other[(this.start - 1):(this.start + nlvls - 3)] <- diff(c(0, start.this)[order(c(0, start.this))])
						attr(allref[[term]], "type") <- 3
					}					
				} else {
					allref[[term]][[1]] <- lvls
					start.new.other[(this.start - 1):(this.start + nlvls - 3)] <- diff(c(0, start.this))
					attr(allref[[term]], "type") <- 3
				}
				this.start <- this.start + nlvls - 1
			} else if (family$family != "binomial" & !monotonic[term]) {
				allref[[term]] <- as.list(lvls)
				attr(allref[[term]], "type") <- 2
			} else {
				if (monotonic[term]) allref[[term]][[1]] <- lvls
				else allref[[term]] <- combinat::permn(lvls)
				attr(allref[[term]], "type") <- 3
			}
		}
    }
    list(allref = allref, terms = t, data = data, monotonic = monotonic, start.new = c(start.new.int, start.new.other, start.new.scale))
}