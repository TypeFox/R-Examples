betadisperF <-
function (d, group, type = c("median", "centroid"), bias.adjust = FALSE) 
{
	#modification of betadisper function (Gavin L. Simpson; bias correction by Adrian Stier and Ben Bolker.) from the Package vegan v2.0-4 (Oksanen & Simpson)


ordimedian<-function (ord, groups, display = "sites", label = FALSE, ...) 
{
    medfun <- function(x, ord) sum(sqrt(rowSums(sweep(ord, 2, 
        x)^2)), na.rm = TRUE)
    dmedfun <- function(x, ord) {
        up <- -sweep(ord, 2, x)
        dn <- sqrt(rowSums(sweep(ord, 2, x)^2))
        colSums(sweep(up, 1, dn, "/"))
    }
    pts <- scores(ord, display = display, ...)
    inds <- names(table(groups))
    medians <- matrix(NA, nrow = length(inds), ncol = ncol(pts))
    rownames(medians) <- inds
    colnames(medians) <- colnames(pts)
    for (i in inds) {
        X <- pts[groups == i, , drop = FALSE]
        if (NROW(X) > 0) 
            medians[i, ] <- optim(apply(X, 2, median, na.rm = TRUE), 
                fn = medfun, gr = dmedfun, ord = X, method = "BFGS")$par
        if (label) 
            ordiArgAbsorber(medians[i, 1], medians[i, 2], label = i, 
                FUN = text, ...)
    }
    invisible(medians)
}

	
ordiArgAbsorber<-function(..., shrink, origin, scaling, triangular, display, 
    choices, const, FUN) 
match.fun(FUN)(...)	

"scores" <-function(x, ...) UseMethod("scores")

	
"scores.default" <-
    function (x, choices, display = c("sites", "species"), ...) 
{
    display <- match.arg(display)
    att <- names(x)
    if (is.data.frame(x) && all(sapply(x, is.numeric)))
        x <- as.matrix(x)
    if (is.list(x) && display == "sites") {
        if ("points" %in% att) 
            X <- x$points
        else if ("rproj" %in% att) 
            X <- x$rproj
        else if ("x" %in% att) 
            X <- x$x
        else if ("scores" %in% att) 
            X <- x$scores
        else if ("sites" %in% att)
            X <- x$sites
        else if("li" %in% att)
            X <- x$li
        else if("l1" %in% att)
            X <- x$l1
        else stop("Can't find scores")
    }
    else if (is.list(x) && display == "species") {
        if ("species" %in% att)
            X <- x$species
        else if ("cproj" %in% att) 
            X <- x$cproj
        else if ("rotation" %in% att) 
            X <- x$rotation
        else if ("loadings" %in% att) 
            X <- x$loadings
        else if ("co" %in% att)
            X <- x$co
        else if ("c1" %in% att)
            X <- x$c1
        else stop("Can't find scores")
    }
    else if (is.numeric(x)) {
        X <- as.matrix(x)
     }
    if (is.null(rownames(X))) {
        root <- substr(display, 1, 4)
        rownames(X) <- paste(root, 1:nrow(X), sep = "")
    }
    if (is.null(colnames(X))) 
        colnames(X) <- paste("Dim", 1:ncol(X), sep = "")
    if (!missing(choices)) {
        choices <- choices[choices <= ncol(X)]
        X <- X[, choices, drop = FALSE]
    }
    X <- as.matrix(X)
    X
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	    dblcen <- function(x, na.rm = TRUE) {
        cnt <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2L, cnt, check.margin = FALSE)
        cnt <- rowMeans(x, na.rm = na.rm)
        sweep(x, 1L, cnt, check.margin = FALSE)
    }


    
    spatialMed <- function(vectors, group, pos) {
        axes <- seq_len(NCOL(vectors))
        spMedPos <- ordimedian(vectors, group, choices = axes[pos])
        spMedNeg <- ordimedian(vectors, group, choices = axes[!pos])
        cbind(spMedPos, spMedNeg)
    }
    Resids <- function(x, c) {
        if (is.matrix(c)) 
            d <- x - c
        else d <- sweep(x, 2, c)
        rowSums(d^2)
    }
    TOL <- 1e-07
    if (!inherits(d, "dist")) 
        stop("distances 'd' must be a 'dist' object")
    if (any(d < -TOL, na.rm = TRUE)) 
        stop("dissimilarities 'd' must be non-negative")
    if (missing(type)) 
        type <- "median"
    type <- match.arg(type)
    if (!is.factor(group)) {
        group <- as.factor(group)
    }
    else {
        group <- droplevels(group)
    }
    n <- attr(d, "Size")
    x <- matrix(0, ncol = n, nrow = n)
    x[row(x) > col(x)] <- d^2
    labs <- attr(d, "Labels")
    if (any(gr.na <- is.na(group))) {
        group <- group[!gr.na]
        x <- x[!gr.na, !gr.na]
        n <- n - sum(gr.na)
        labs <- labs[!gr.na]
        warning("Missing observations due to 'group' removed.")
    }
    if (any(x.na <- apply(x, 1, function(x) any(is.na(x))))) {
        x <- x[!x.na, !x.na]
        group <- group[!x.na]
        n <- n - sum(x.na)
        labs <- labs[!x.na]
        warning("Missing observations due to 'd' removed.")
    }
    x <- x + t(x)
 x<- dblcen(x)
    e <- eigen(-x/2, symmetric = TRUE)
    vectors <- e$vectors
    eig <- e$values
    eig <- eig[(want <- abs(eig/eig[1]) > TOL)]
    vectors <- vectors[, want, drop = FALSE] %*% diag(sqrt(abs(eig)), 
        nrow = length(eig))
    pos <- eig > 0
    centroids <- switch(type, centroid = apply(vectors, 2, function(x) tapply(x,group, mean)), median = spatialMed(vectors, group, pos))
    if(is.matrix(centroids)==FALSE){
    centroids<-t(as.matrix(centroids))
    rownames(centroids)<-"YES"}
    dist.pos <- Resids(vectors[, pos, drop = FALSE], centroids[group, 
        pos, drop = FALSE])
    dist.neg <- 0
    if (any(!pos)) 
        dist.neg <- Resids(vectors[, !pos, drop = FALSE], centroids[group, 
            !pos, drop = FALSE])
    zij <- sqrt(abs(dist.pos - dist.neg))
    if (bias.adjust) {
        n.group <- table(group)
        zij <- zij * sqrt(n.group[group]/(n.group[group] - 1))
    }
    colnames(vectors) <- names(eig) <- paste("PCoA", seq_along(eig), 
        sep = "")
    if (is.matrix(centroids)) 
        colnames(centroids) <- names(eig)
    else names(centroids) <- names(eig)
    rownames(vectors) <- names(zij) <- labs
    retval <- list(eig = eig, vectors = vectors, distances = zij, 
        group = group, centroids = centroids, call = match.call())
    class(retval) <- "betadisper"
    attr(retval, "method") <- attr(d, "method")
    attr(retval, "type") <- type
    attr(retval, "bias.adjust") <- bias.adjust
    retval
}
