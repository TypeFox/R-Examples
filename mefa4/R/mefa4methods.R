## accessor methods 

## for new mefa classes
setGeneric("xtab", function(x) standardGeneric("xtab"))
setGeneric("samp", function(x) standardGeneric("samp"))
setGeneric("taxa", function(x) standardGeneric("taxa"))
setMethod("xtab", signature(x = "Mefa"), function(x) x@xtab)
setMethod("samp", signature(x = "Mefa"), function(x) x@samp)
setMethod("taxa", signature(x = "Mefa"), function(x) x@taxa)

## for old mefa classes
setMethod("xtab", signature(x = "mefa"), function(x) x$xtab)
setMethod("samp", signature(x = "mefa"), function(x) x$samp)
setMethod("taxa", signature(x = "mefa"), function(x) x$taxa)

## setters and replacement

setGeneric("xtab<-", function(x, value) standardGeneric("xtab<-"))
setGeneric("samp<-", function(x, value) standardGeneric("samp<-"))
setGeneric("taxa<-", function(x, value) standardGeneric("taxa<-"))

setReplaceMethod("xtab", signature(x = "Mefa", value = "MefaMatrix"),
    function(x, value) {
        value <- as(value, "dgCMatrix")
        if (x@join == "left") {
            rkeep <- rownames(value)
            ckeep <- colnames(value)
            x@xtab <- value
            if (!is.null(x@samp))
                x@samp <- x@samp[rkeep,,drop=FALSE]
            if (!is.null(x@taxa))
                x@taxa <- x@taxa[ckeep,,drop=FALSE]
        }
        if (x@join == "inner") {
            rkeep <- if (!is.null(x@samp)) {
                intersect(rownames(value), rownames(x@samp))
            } else {
                rownames(value)
            }
            ckeep <- if (!is.null(x@taxa)) {
                intersect(colnames(value), rownames(x@taxa))
            } else {
                colnames(value)
            }
#            XTAB <- value[rkeep, ckeep,drop=FALSE]
#            if (is.null(dim(XTAB))) {
#                dim(XTAB) <- c(length(rkeep), length(ckeep))
#                dimnames(XTAB) <- list(rkeep, ckeep)
#            }
#            x@xtab <- XTAB
            x@xtab <- value[rkeep, ckeep,drop=FALSE]
            if (!is.null(x@samp))
                x@samp <- x@samp[rkeep,,drop=FALSE]
            if (!is.null(x@taxa))
                x@taxa <- x@taxa[ckeep,,drop=FALSE]
        }
        if (!is.null(x@samp))
            rownames(x@samp) <- rkeep
        if (!is.null(x@taxa))
            rownames(x@taxa) <- ckeep
        x
})
setReplaceMethod("samp", signature(x = "Mefa", value = "MefaDataFrame"), 
    function(x, value) {
        if (!is.null(value)) {
            if (x@join == "left") {
                rkeep <- rownames(x@xtab)
                x@samp <- value[rkeep,,drop=FALSE]
            }
            if (x@join == "inner") {
                rkeep <- intersect(rownames(x@xtab), rownames(value))
                x@xtab <- x@xtab[rkeep,,drop=FALSE]
                x@samp <- value[rkeep,,drop=FALSE]
            }
            rownames(x@samp) <- rkeep
        } else x@samp <- NULL
        x
})
setReplaceMethod("taxa", signature(x = "Mefa", value = "MefaDataFrame"), 
    function(x, value) {
        if (!is.null(value)) {
            if (x@join == "left") {
                ckeep <- colnames(x@xtab)
                x@taxa <- value[ckeep,,drop=FALSE]
            }
            if (x@join == "inner") {
                ckeep <- intersect(colnames(x@xtab), rownames(value))
                x@xtab <- x@xtab[,ckeep,drop=FALSE]
                x@taxa <- value[ckeep,,drop=FALSE]
            }
            rownames(x@taxa) <- ckeep
        } else x@taxa <- NULL
        x
})

## subsetting [
## TODO: vary for different signatures???

setMethod("[", signature(x = "Mefa", i = "ANY", 
        j = "ANY", drop = "ANY"),
    function(x, i, j, ..., drop) {
        if (missing(i))
            i <- 1:dim(x)[1]
        if (missing(j))
            j <- 1:dim(x)[2]
        if (missing(drop))
            drop <- FALSE
        if (any(is.na(i)))
            stop("index contains 'NA'")
        if (any(is.na(j)))
            stop("index contains 'NA'")
#        XTAB <- x@xtab[i,j,drop=FALSE]
#        if (is.null(dim(XTAB))) {
#            li <- if (is.logical(i))
#                sum(i) else length(i)
#            lj <- if (is.logical(j))
#                sum(j) else length(j)
#            dim(XTAB) <- c(li, lj)
#            dimnames(XTAB) <- list(rownames(x@xtab)[i], colnames(x@xtab)[j])
#        }
#        x@xtab <- XTAB
        x@xtab <- x@xtab[i,j,drop=FALSE]
        if (!is.null(x@samp)) {
            x@samp <- x@samp[i,,drop=FALSE]
            if (drop)
                x@samp <- lapply(x@samp, function(z) z[drop=TRUE])
        }
        if (!is.null(x@taxa)) {
            x@taxa <- x@taxa[j,,drop=FALSE]
            if (drop)
                x@taxa <- lapply(x@taxa, function(z) z[drop=TRUE])
        }
        x
})

## coercion

setMethod("as.matrix", "Mefa", function(x) as.matrix(x@xtab))
setAs(from = "matrix", to = "Mefa", def = function(from) Mefa(from))
setAs(from = "Mefa", to = "sparseMatrix", def = function(from) from@xtab)
setAs(from = "sparseMatrix", to = "Mefa", def = function(from) Mefa(from))

## general methods

setMethod("dim", "Mefa", function(x) dim(x@xtab))
setMethod("dimnames", "Mefa", function(x) dimnames(x@xtab))

setMethod("dimnames<-", signature(x = "Mefa", value = "list"), 
    function(x, value) {
        dimnames(x@xtab) <- value
        if (!is.null(x@samp))
            rownames(x@samp) <- value[[1]]
        if (!is.null(x@taxa))
            rownames(x@taxa) <- value[[2]]
        x
})

## transpose, why not?
setMethod("t", "Mefa", function(x) {
    new("Mefa", xtab = t(x@xtab),
        samp = x@taxa, taxa = x@samp,
        join = x@join)
})

## show for Mefa
setMethod("show", "Mefa", function(object) {
    d <- dim(object)
    cat("Object of class \"Mefa\"\n")
    cat("  ..@ xtab:", d[1], "x", d[2], "sparse Matrix\n")
    if (!is.null(object@samp)) {
        cat("  ..@ samp: data frame with", ncol(object@samp), "variables\n")
    } else {
        cat("  ..@ samp: NULL\n")
    }
    if (!is.null(object@taxa)) {
        cat("  ..@ taxa: data frame with", ncol(object@taxa), "variables\n")
    } else {
        cat("  ..@ taxa: NULL\n")
    }
    cat("  ..@ join:", object@join, "\n")
    invisible(object)
})

setMethod("stack", "Mefa", function(x, ...) {
    d <- dim(x)
    dn <- dimnames(x)
    X <- data.frame(
        samp=as.factor(rep(dn[[1]], d[2])),
        taxa=as.factor(rep(dn[[2]], each=d[1])),
        values=as.numeric(xtab(x)))
    SAMP <- samp(x)
    if (!is.null(SAMP)) {
        colnames(SAMP) <- paste("samp", colnames(SAMP), sep="_")
        X <- data.frame(X, 
            SAMP[match(X$samp, rownames(SAMP)),])
    }
    TAXA <- taxa(x)
    if (!is.null(TAXA)) {
        colnames(TAXA) <- paste("taxa", colnames(TAXA), sep="_")
        X <- data.frame(X, 
            TAXA[match(X$taxa, rownames(TAXA)),])
    }
    X
})
