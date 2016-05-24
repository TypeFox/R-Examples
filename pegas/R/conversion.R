## conversion.R (2016-03-16)

##   Conversion Among Allelic Data Classes

## Copyright 2009-2016 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

loci2genind <- function(x)
{
    ipop <- which(names(x) == "population")
    pop <- if (length(ipop)) x[, ipop] else NULL
    df2genind(as.matrix(x[, attr(x, "locicol")]), sep = "[/\\|]",
              pop = pop, NA.char = ".") # fix by Thibaut (2015-11-10)
}

as.loci <- function(x, ...) UseMethod("as.loci")

as.loci.genind <- function(x, ...)
{
    obj <- genind2df(x, sep = "/")
    icol <- 1:ncol(obj)
    pop <- which(names(obj) == "pop")
    if (length(pop)) {
        names(obj)[pop] <- "population"
        icol <- icol[-pop]
    }
    for (i in icol) obj[, i] <- factor(obj[, i] )
    class(obj) <- c("loci", "data.frame")
    attr(obj, "locicol") <- icol
    .check.order.alleles(obj)
}

genind2loci <- function(x) as.loci.genind(x)

## to be sure that alleles are sorted in their ASCII code
## (not in lexicographical order) whatever the locale
.sort.alleles <- function(x, index.only = FALSE)
{
    locale <- Sys.getlocale("LC_COLLATE")
    if (!identical(locale, "C")) {
        Sys.setlocale("LC_COLLATE", "C")
        on.exit(Sys.setlocale("LC_COLLATE", locale))
    }
    if (index.only) order(x) else sort(x)
}

.check.order.alleles <- function(x)
{
    reorder.alleles <- function(x) {
        for (i in seq_along(x)) {
            y <- x[i]
            if (!length(grep("/", y))) next # phased genotype (mixed with unphased ones)
            y <- unlist(strsplit(y, "/"))
            y <- paste(.sort.alleles(y), collapse = "/")
            x[i] <- y
        }
        x
    }

    for (k in attr(x, "locicol")) {
        y <- x[, k]
        if (is.numeric(y)) { # haploid with alleles coded with numerics
            x[, k] <- factor(y)
            next
        }
        lv <- levels(y)
        if (!length(grep("/", lv))) next # if haploid or phased genotype
        a <- reorder.alleles(lv) # works with all levels of ploidy > 1
        if (!identical(a, lv)) levels(x[, k]) <- a
    }
    x
}

as.loci.data.frame <-
    function(x, allele.sep = "/|", col.pop = NULL, col.loci = NULL, ...)
{
    if (is.null(col.pop)) {
        ipop <- which(tolower(names(x)) == "population")
        if (length(ipop)) col.pop <- ipop
    }
    if (is.character(col.pop))
        col.pop <- which(names(x) == col.pop)
    if (is.numeric(col.pop)) {
        names(x)[col.pop] <- "population"
        x[, col.pop] <- factor(x[, col.pop])
    }
    if (is.null(col.loci)) {
        col.loci <- 1:ncol(x)
        if (is.numeric(col.pop))
            col.loci <- col.loci[-col.pop]
    }
    if (is.character(col.loci))
        col.loci <- match(col.loci, names(x))
    if (allele.sep != "/|") {
        if (allele.sep == "")
            stop("alleles within a genotype must be separated")
        for (i in col.loci)
            levels(x[, i]) <- gsub(allele.sep, "/", levels(x[, i]))
    }
    class(x) <- c("loci", "data.frame")
    attr(x, "locicol") <- col.loci
    .check.order.alleles(x)
}

as.loci.factor <- function(x, allele.sep = "/|", ...)
    as.loci.data.frame(data.frame(x), allele.sep = allele.sep, ...)

as.loci.character <- function(x, allele.sep = "/|", ...)
    as.loci.data.frame(data.frame(factor(x)), allele.sep = allele.sep, ...)

alleles2loci <- function(x, ploidy = 2, rownames = NULL, population = NULL,
                         phased = FALSE)
{
    withPop <- !is.null(population)
    x <- as.data.frame(x)
    if (is.null(rownames)) {
        idx <- rownames(x)
        if (is.null(idx)) idx <- as.character(seq_len(d[1]))
    } else {
        idx <- as.character(x[[rownames]])
        x[[rownames]] <- NULL
        if (withPop && rownames < population)
            population <- population - 1
    }
    if (withPop) {
        pop <- x[, population]
        x <- x[, -population]
    }
    d <- dim(x)
    if (d[2] %% ploidy) stop("number of columns not a multiple of ploidy")
    nloci <- d[2] / ploidy
    start <- seq(1, by = ploidy, length.out = nloci)
    end <- start + ploidy - 1
    loci.nms <- colnames(x)[start]
    obj <- vector("list", nloci)
    sep <- if (phased) "|" else "/"
    foo <- function(...) paste(..., sep = sep)
    for (i in seq_len(nloci))
        obj[[i]] <- factor(do.call(foo, x[, start[i]:end[i]]))
    names(obj) <- loci.nms
    obj <- as.data.frame(obj, row.names = idx)
    obj <- as.loci(obj)
    if (withPop) obj$population <- factor(pop)
    obj
}

na.omit.loci <- function(object, na.alleles = c("0", "."), ...)
{
    pat <- c(paste0("^", na.alleles, "/"), paste0("/", na.alleles, "$"), paste0("/", na.alleles, "/"))
    pat <- paste(pat, collapse = "|")
    drop <- logical(nrow(object))
    M <- 1:ncol(object)
    for (i in attr(object, "locicol")) {
        x <- object[[i]]
        if (length(na <- grep(pat, x))) drop[na] <- TRUE
        if (any(na <- is.na(x))) drop[na] <- TRUE
    }
    object <- object[!drop, ]
    for (i in M) {
        if (is.factor(x <- object[[i]])) {
            drop <- tabulate(x, nlevels(x)) == 0
            if (any(drop)) object[[i]] <- factor(x)
        }
    }
    object
}
