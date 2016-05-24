## TODO: mbind for matrix could use sparse matrix code
## check for necessity of [,,drop=FALSE] in mbind for 1x cases
## check problem with merge for samp/taxa

## groupSums and groupMeans

setGeneric("groupSums", function(object, ...) standardGeneric("groupSums"))

## MARGIN indicates what to group (1: group rows, 2: group cols)
setMethod("groupSums", "matrix", function(object, MARGIN, by, na.rm = FALSE, ...) {
    if (any(is.na(by)))
        stop("'NA' not allowed in 'by'")
    if (any(is.na(object)) && !na.rm)
        stop("'NA' found in 'object'")
    if (!(MARGIN %in% 1:2))
        stop("'MARGIN' must be in 1:2")
    if (length(MARGIN) != 1)
        stop("MARGIN = 1:2 not yet implemented")
    if (length(by) != dim(object)[MARGIN])
        stop("Non conforming 'object', 'MARGIN' and 'by'")
    mm <- as(factor(by, levels=unique(by)), "sparseMatrix")
    rownames(mm) <- unique(by)
    object <- as(object, "sparseMatrix")
    if (na.rm)
        object[is.na(object)] <- 0
    if (MARGIN == 2) {
        out <- t(mm %*% t(object))
    } else {
        out <- mm %*% object
    }
    as.matrix(out)
})
setMethod("groupSums", "sparseMatrix", function(object, MARGIN, by, na.rm = FALSE, ...) {
    if (any(is.na(by)))
        stop("'NA' not allowed in 'by'")
    if (any(is.na(object)) && !na.rm)
        stop("'NA' found in 'object'")
    if (!(MARGIN %in% 1:2))
        stop("'MARGIN' must be in 1:2")
    if (length(MARGIN) != 1)
        stop("MARGIN = 1:2 not yet implemented")
    if (length(by) != dim(object)[MARGIN])
        stop("Non conforming 'object', 'MARGIN' and 'by'")
    mm <- as(factor(by, levels=unique(by)), "sparseMatrix")
    rownames(mm) <- unique(by)
    if (na.rm)
        object[is.na(object)] <- 0
    if (MARGIN == 2) {
        out <- t(mm %*% t(object))
    } else {
        out <- mm %*% object
    }
    out
})
## replace is a replacement object for the affected non xtab slot (samp, taxa)
setMethod("groupSums", "Mefa", function(object, MARGIN, by, replace, na.rm = FALSE, ...) {
    x <- groupSums(object@xtab, MARGIN, by, na.rm, ...)
    if (missing(replace))
        replace <- NULL
    JOIN <- object@join
    if (MARGIN == 2) {
#        if (is.null(object@samp) && is.null(replace))
#            JOIN <- "inner"
        new("Mefa", xtab = x, samp = object@samp,
            taxa = replace, join = JOIN)
    } else {
#        if (is.null(object@taxa) && is.null(replace))
#            JOIN <- "inner"
        new("Mefa", xtab = x, samp = replace,
            taxa = object@taxa, join = JOIN)
    }
})

setGeneric("groupMeans", function(object, ...) standardGeneric("groupMeans"))

setMethod("groupMeans", "sparseMatrix", function(object, MARGIN, by, na.rm = FALSE, ...) {
    x <- groupSums(object, MARGIN, by, na.rm, ...)
    if (any(is.na(x))) {
        n <- groupSums(as(!is.na(object), "sparseMatrix"), MARGIN, by, na.rm, ...)
        out <- x / n
    } else {
        n <- rowSums(as(factor(by, levels=unique(by)), "sparseMatrix"))
        out <- if (MARGIN == 1)
            x / n else t(t(x) / n)
    }
    as(out, "sparseMatrix")
})
setMethod("groupMeans", "matrix", function(object, MARGIN, by, na.rm = FALSE, ...) {
    as.matrix(groupMeans(as(object, "sparseMatrix"), MARGIN, by, na.rm, ...))
})
setMethod("groupMeans", "Mefa", function(object, MARGIN, by, replace, na.rm = FALSE, ...) {
    x <- groupMeans(as(object, "sparseMatrix"), MARGIN, by, na.rm, ...)
    if (missing(replace))
        replace <- NULL
    JOIN <- object@join
    if (MARGIN == 2) {
#        if (is.null(object@samp) && is.null(replace))
#            JOIN <- "inner"
        new("Mefa", xtab = x, samp = object@samp,
            taxa = replace, join = JOIN)
    } else {
#        if (is.null(object@taxa) && is.null(replace))
#            JOIN <- "inner"
        new("Mefa", xtab = x, samp = replace,
            taxa = object@taxa, join = JOIN)
    }
})

## mbind: joining 2 Mefa objects

setGeneric("mbind", function(x, y, fill, ...) standardGeneric("mbind"))

setMethod("mbind", signature(x="matrix", y="matrix", fill="ANY"), 
    function(x, y, fill, ...) {
        if (missing(fill))
            fill <- NA
        if (length(x) == 0)
            stop("length of 'x' must not be 0")
        if (length(y) == 0)
            stop("length of 'y' must not be 0")
        if (is.null(dimnames(x)))
            stop("dimnames of 'x' must not be NULL")
        if (is.null(dimnames(y)))
            stop("dimnames of 'y' must not be NULL")
        r1 <- rownames(x)
        c1 <- colnames(x)
        r2 <- rownames(y)
        c2 <- colnames(y)
        if (setequal(r1, r2) && setequal(c1, c2))
            return(x)
        rx <- setdiff(r1, r2)
        rxy <- intersect(r1, r2)
        ry <- setdiff(r2, r1)
        cx <- setdiff(c1, c2)
        cxy <- intersect(c1, c2)
        cy <- setdiff(c2, c1)
        xx <- x[c(rx, rxy), c(cx, cxy), drop=FALSE]
        z1 <- matrix(fill, length(ry), length(cx))
        z2 <- matrix(fill, length(rx), length(cy))
        yx1 <- y[ry, cxy, drop=FALSE]
        yx2 <- y[rxy, cy, drop=FALSE]
        yy <- y[ry, cy, drop=FALSE]
        if (length(ry) > 0) {
            part1 <- cbind(z1, yx1)
            part2 <- rbind(xx, part1)
        } else part2 <- xx
        if (length(cy) > 0) {
            part3 <- rbind(z2, yx2)
            part4 <- rbind(part3, yy)
            part5 <- cbind(part2, part4)
        } else part5 <- part2
        rownames(part5) <- c(rx, rxy, ry)
        colnames(part5) <- c(cx, cxy, cy)
        part5
})
setMethod("mbind", signature(x="sparseMatrix", y="sparseMatrix", fill="ANY"), 
    function(x, y, fill, ...) {
        if (missing(fill))
            fill <- NA
        if (length(x) == 0)
            stop("length of 'x' must not be 0")
        if (length(y) == 0)
            stop("length of 'y' must not be 0")
        if (is.null(dimnames(x)))
            stop("dimnames of 'x' must not be NULL")
        if (is.null(dimnames(y)))
            stop("dimnames of 'y' must not be NULL")
        r1 <- rownames(x)
        c1 <- colnames(x)
        r2 <- rownames(y)
        c2 <- colnames(y)
        if (setequal(r1, r2) && setequal(c1, c2))
            return(x)
        rx <- setdiff(r1, r2)
        rxy <- intersect(r1, r2)
        ry <- setdiff(r2, r1)
        cx <- setdiff(c1, c2)
        cxy <- intersect(c1, c2)
        cy <- setdiff(c2, c1)
        xx <- x[c(rx, rxy), c(cx, cxy), drop=FALSE]
        z1 <- as(matrix(fill, length(ry), length(cx)), "sparseMatrix")
        z2 <- as(matrix(fill, length(rx), length(cy)), "sparseMatrix")
        yx1 <- y[ry, cxy, drop=FALSE]
        yx2 <- y[rxy, cy, drop=FALSE]
        yy <- y[ry, cy, drop=FALSE]
        if (length(ry) > 0) {
            part1 <- cBind(z1, yx1)
            part2 <- rBind(xx, part1)
        } else part2 <- xx
        if (length(cy) > 0) {
            part3 <- rBind(z2, yx2)
            part4 <- rBind(part3, yy)
            part5 <- cBind(part2, part4)
        } else part5 <- part2
        rownames(part5) <- c(rx, rxy, ry)
        colnames(part5) <- c(cx, cxy, cy)
        part5
})
setMethod("mbind", signature(x="Mefa", y="Mefa", fill="ANY"), 
    function(x, y, fill, drop, ...) {
        if (missing(drop))
            drop <- FALSE
        if (missing(fill))
            fill <- NA
        ## xtab
        xtabx <- x@xtab
        xtaby <- y@xtab
        if (length(xtabx) == 0)
            stop("length of 'x' must not be 0")
        if (length(xtaby) == 0)
            stop("length of 'y' must not be 0")
        r1 <- rownames(xtabx)
        c1 <- colnames(xtabx)
        r2 <- rownames(xtaby)
        c2 <- colnames(xtaby)
        if (setequal(r1, r2) && setequal(c1, c2))
            return(x)
        rx <- setdiff(r1, r2)
        rxy <- intersect(r1, r2)
        ry <- setdiff(r2, r1)
        cx <- setdiff(c1, c2)
        cxy <- intersect(c1, c2)
        cy <- setdiff(c2, c1)
        xx <- xtabx[c(rx, rxy), c(cx, cxy), drop=FALSE]
        z1 <- as(matrix(fill, length(ry), length(cx)), "sparseMatrix")
        z2 <- as(matrix(fill, length(rx), length(cy)), "sparseMatrix")
        yx1 <- xtaby[ry, cxy, drop=FALSE]
        yx2 <- xtaby[rxy, cy, drop=FALSE]
        yy <- xtaby[ry, cy, drop=FALSE]
        if (length(ry) > 0) {
            part1 <- cBind(z1, yx1)
            part2 <- rBind(xx, part1)
        } else part2 <- xx
        if (length(cy) > 0) {
            part3 <- rBind(z2, yx2)
            part4 <- rBind(part3, yy)
            part5 <- cBind(part2, part4)
        } else part5 <- part2
        rownames(part5) <- c(rx, rxy, ry)
        colnames(part5) <- c(cx, cxy, cy)
        ## samp
        sampx <- x@samp
        sampy <- y@samp
        if (is.null(sampx) && is.null(sampy))
            sm2 <- NULL
        if (!is.null(sampx) && !is.null(sampy)) {
            sampx2 <- data.frame(ROWNAMES=rownames(sampx), SAMPPART=1, sampx[,,drop=FALSE])
            sampy2 <- data.frame(ROWNAMES=rownames(sampy), SAMPPART=2, sampy[,,drop=FALSE])
            saa <- intersect(names(sampy2), names(sampx2))
            sbb <- c(saa, setdiff(names(sampx2), names(sampy2)))
            sm <- merge(sampx2[,sbb, drop=FALSE], sampy2[,saa, drop=FALSE], all=TRUE)
            sid1 <- which(sm$SAMPPART==1 & sm$ROWNAMES %in% r1)
            sid2 <- which(sm$SAMPPART==2 & sm$ROWNAMES %in% ry)
            sm2 <- sm[c(sid1, sid2), , drop=FALSE]
            rownames(sm2) <- sm2$ROWNAMES
            scc <- setdiff(names(sampy), names(sampx))
            sm2 <- data.frame(sm2, sampy[rownames(sm2),scc,drop=FALSE])
            sm2$ROWNAMES <- sm2$SAMPPART <- NULL
        }
        if (!is.null(sampx) && is.null(sampy))
            sm2 <- sampx
        if (is.null(sampx) && !is.null(sampy))
            sm2 <- sampy[c(r1,ry), , drop=FALSE]
        ## taxa
        taxax <- x@taxa
        taxay <- y@taxa
        if (is.null(taxax) && is.null(taxay))
            tm2 <- NULL
        if (!is.null(taxax) && !is.null(taxay)) {
            taxax2 <- data.frame(ROWNAMES=rownames(taxax), TAXAPART=1, taxax[,,drop=FALSE])
            taxay2 <- data.frame(ROWNAMES=rownames(taxay), TAXAPART=2, taxay[,,drop=FALSE])
            taa <- intersect(names(taxay2), names(taxax2))
            tbb <- c(taa, setdiff(names(taxax2), names(taxay2)))
            tm <- merge(taxax2[,tbb, drop=FALSE], taxay2[,taa, drop=FALSE], all=TRUE)
            tid1 <- which(tm$TAXAPART==1 & tm$ROWNAMES %in% c1)
            tid2 <- which(tm$TAXAPART==2 & tm$ROWNAMES %in% cy)
            tm2 <- tm[c(tid1, tid2), , drop=FALSE]
            rownames(tm2) <- tm2$ROWNAMES
            tcc <- setdiff(names(taxay), names(taxax))
            tm2 <- data.frame(tm2, taxay[rownames(tm2),tcc,drop=FALSE])
            tm2$ROWNAMES <- tm2$TAXAPART <- NULL
        }
        if (!is.null(taxax) && is.null(taxay))
            tm2 <- taxax
        if (is.null(taxax) && !is.null(taxay))
            tm2 <- taxay[c(c1,cy), , drop=FALSE]
        ## assembling
        Mefa(part5, sm2, tm2, join="left", drop)
})

setGeneric("mbind2", function(x, y, fill, ...) standardGeneric("mbind2"))

setMethod("mbind2", signature(x="matrix", y="matrix", fill="ANY"), 
    function(x, y, fill, ...) {
        as.matrix(mbind2(x=as(x, "sparseMatrix"),
            y=as(y, "sparseMatrix"),
            fill=fill, ...))
})

setMethod("mbind2", signature(x="sparseMatrix", y="sparseMatrix", fill="ANY"), 
    function(x, y, fill, ...) {
        if (missing(fill))
            fill <- NA
        if (length(x) == 0)
            stop("length of 'x' must not be 0")
        if (length(y) == 0)
            stop("length of 'y' must not be 0")
        if (is.null(dimnames(x)))
            stop("dimnames of 'x' must not be NULL")
        if (is.null(dimnames(y)))
            stop("dimnames of 'y' must not be NULL")
        rn <- union(rownames(x), rownames(y))
        cn <- union(colnames(x), colnames(y))
        rn0 <- intersect(rownames(x), rownames(y))
        cn0 <- intersect(colnames(x), colnames(y))
        out <- Matrix(fill, 
            length(rn), length(cn), 
            dimnames=list(rn, cn),
            sparse=TRUE)
        out <- as(out, "dgCMatrix")
        out[rownames(x), colnames(x)] <- x
        out[rownames(y), colnames(y)] <- y
        out[rn0, cn0] <- x[rn0, cn0] + y[rn0, cn0]
        out
})
