# Copyright 2004-5 (c) Roger Bivand

dotsInPolys <- function(pl, x, f="random", offset, compatible=FALSE) {
    if (!is.character(f)) stop("f must be a character string")
    if (f != "random" && f != "regular") stop(paste(f, "not supported"))
#    if (inherits(pl, "polylist")) pl <- .polylist2SpP(pl)
    if (!is(pl, "SpatialPolygons")) stop("unknown class of input polygons")
    pls <- slot(pl, "polygons")
    IDs <- sapply(pls, function(i) slot(i, "ID"))
    if (length(pls) != length(x)) stop("different lengths")
    if (!inherits(x, "integer")) {
        x <- as.integer(x)
        warning("x coerced to integer")
    }
    n <- length(pls)
    if (n < 1) stop("zero Polygons")
    res <- vector(mode="list", length=n)
    ID_out <- NULL
    if (missing(offset)) {
	if (f == "random") offset <- runif(2)
	else offset <- c(0.5,0.5)
    }
    for (i in 1:n) {
        if (x[i] > 0) {
		#EJP: deprecate sample.Polygons
	    #ires <- sample.Polygons(pls[[i]], x[i], type=f, offset=offset)
	    ires <- spsample(pls[[i]], x[i], type=f, offset=offset)
	    if (!is.null(ires)) res[[i]] <- ires
	    if (!is.null(res[[i]])) ID_out <- c(ID_out, IDs[i])
	}
    }
    if (!compatible) {
        resa <- do.call("rbind", lapply(res, function(x) 
	    if (!is.null(x)) coordinates(x)))
	reps <- unlist(sapply(res, function(x) 
	    if (!is.null(x)) nrow(coordinates(x))))
	res <- data.frame(resa, rep(ID_out, reps))
	names(res) <- c("x", "y", "ID")
	coordinates(res) <- c("x", "y")
    } else {
	j <- 1
        for (i in 1:n) {
	    if (!is.null(res[[i]])) {
		res[[i]] <- coordinates(res[[i]])
		attr(res[[i]], "ID") <- ID_out[j]
		j <- j+1
	    }
        }
    }
    res
}

symbolsInPolys <- function(pl, dens, symb="+", compatible=FALSE) {
#    if (inherits(pl, "polylist")) pl <- .polylist2SpP(pl)
    if (!is(pl, "SpatialPolygons")) stop("unknown class of input polygons")
    if (!is(pl, "SpatialPolygons")) stop("unknown class of input polygons")
    pls <- slot(pl, "polygons")
    n <- length(pls)
    if (n < 1) stop("zero Polygons")
    if (n != length(dens)) dens <- rep(dens[1], n)
    if (n != length(symb)) symb <- rep(symb[1], n)
    areas <- lapply(pls, function(x) sapply(slot(x, "Polygons"),
        function(i) slot(i, "area")))
    holes <- lapply(pls, function(x) sapply(slot(x, "Polygons"),
        function(i) slot(i, "hole")))
    counts <- vector(mode="list", n)
    for (i in 1:n) {
	cntvec <- NULL
        for (j in 1:length(areas[[i]])) {
	    cntvec[j] <- ifelse(holes[[i]][j], 0, areas[[i]][j] * dens[i])
	}
        counts[[i]] <- as.integer(cntvec)
    }

    res <- vector(mode="list", n)
    symb_out <- NULL
    for (i in 1:n) {
	px <- as.integer(sum(counts[[i]]))
        if (px > 0) {
			#EJP: deprecate sample.Polygons:
            #ires <- sample.Polygons(pls[[i]], px, type="regular", 
            ires <- spsample(pls[[i]], px, type="regular", offset=c(0.5,0.5))
	    if (!is.null(ires)) res[[i]] <- ires
	    if (!is.null(res[[i]])) symb_out <- c(symb_out, symb[i])
	}
    }

    if (!compatible) {
        resa <- do.call("rbind", lapply(res, function(x) 
	    if(!is.null(x)) coordinates(x)))
	reps <- unlist(sapply(res, function(x) 
	    if(!is.null(x)) nrow(coordinates(x))))
	res <- data.frame(resa, rep(symb_out, reps))
	names(res) <- c("x", "y", "symb")
	coordinates(res) <- c("x", "y")
    } else {
	j <- 1
        for (i in 1:n) {
	    if (!is.null(res[[i]])) {
		res[[i]] <- coordinates(res[[i]])
		attr(res[[i]], "symb") <- symb_out[j]
		j <- j+1
	    }
        }
    }
    res
}


