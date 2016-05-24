## DNA.R (2015-11-18)

##   Manipulations and Comparisons of DNA Sequences

## Copyright 2002-2015 Emmanuel Paradis, 2015 Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

DNAbin2indel <- function(x)
{
    if (is.list(x)) x <- as.matrix(x)
    d <- dim(x)
    s <- as.integer(d[2])
    n <- as.integer(d[1])
    res <- .C(DNAbin2indelblock, x, n, s, integer(n*s), NAOK = TRUE)[[4]]
    dim(res) <- d
    rownames(res) <- rownames(x)
    res
}

labels.DNAbin <- function(object, ...)
{
    if (is.list(object)) return(names(object))
    if (is.matrix(object)) return(rownames(object))
    NULL
}

del.gaps <- function(x)
{
    deleteGaps <- function(x) {
        i <- which(x == 4)
        if (length(i)) x[-i] else x
    }

    if (!inherits(x, "DNAbin")) x <- as.DNAbin(x)
    if (is.matrix(x)) {
        n <- dim(x)[1]
        y <- vector("list", n)
        for (i in 1:n) y[[i]] <- x[i, ]
        names(y) <- rownames(x)
        x <- y
        rm(y)
    }
    if (!is.list(x)) return(deleteGaps(x))
    x <- lapply(x, deleteGaps)
    class(x) <- "DNAbin"
    x
}

del.colgapsonly <- function(x)
{
    if (!inherits(x, "DNAbin")) x <- as.DNAbin(x)
    if (!is.matrix(x))
        stop("DNA sequences not in a matrix")
    foo <- function(x) all(x == 4)
    i <- which(apply(x, 2, foo))
    if (length(i)) x <- x[, -i]
    x
}

as.alignment <- function(x)
{
    if (is.list(x)) n <- length(x)
    if (is.matrix(x)) n <- dim(x)[1]
    seq <- character(n)
    if (is.list(x)) {
        nam <- names(x)
        for (i in 1:n)
          seq[i] <- paste(x[[i]], collapse = "")
    }
    if (is.matrix(x)) {
        nam <- dimnames(x)[[1]]
        for (i in 1:n)
          seq[i] <- paste(x[i, ], collapse = "")
    }
    obj <- list(nb = n, seq = seq, nam = nam, com = NA)
    class(obj) <- "alignment"
    obj
}

"[.DNAbin" <- function(x, i, j, drop = FALSE)
{
    ans <- NextMethod("[", drop = drop)
    class(ans) <- "DNAbin"
    ans
}

as.matrix.DNAbin <- function(x, ...)
{
    if (is.matrix(x)) return(x)
    if (is.vector(x)) {
        dim(x) <- c(1, length(x))
        return(x)
    }
    s <- unique(unlist(lapply(x, length)))
    if (length(s) != 1)
        stop("DNA sequences in list not of the same length.")
    nms <- names(x)
    n <- length(x)
    y <- matrix(as.raw(0), n, s)
    for (i in seq_len(n)) y[i, ] <- x[[i]]
    rownames(y) <- nms
    class(y) <- "DNAbin"
    y
}

as.list.DNAbin <- function(x, ...)
{
    if (is.list(x)) return(x)
    if (is.null(dim(x))) obj <- list(x) # cause is.vector() doesn't work
    else { # matrix
        n <- nrow(x)
        obj <- vector("list", n)
        for (i in 1:n) obj[[i]] <- x[i, ]
        names(obj) <- rownames(x)
    }
    class(obj) <- "DNAbin"
    obj
}

rbind.DNAbin <- function(...)
### works only with matrices for the moment
{
    obj <- list(...)
    n <- length(obj)
    if (n == 1) return(obj[[1]])
    for (i in 1:n)
        if (!is.matrix(obj[[1]]))
            stop("the 'rbind' method for \"DNAbin\" accepts only matrices")
    NC <- unlist(lapply(obj, ncol))
    if (length(unique(NC)) > 1)
        stop("matrices do not have the same number of columns.")
    for (i in 1:n) class(obj[[i]]) <- NULL
    for (i in 2:n) obj[[1]] <- rbind(obj[[1]], obj[[i]])
    structure(obj[[1]], class = "DNAbin")
}

cbind.DNAbin <-
    function(..., check.names = TRUE, fill.with.gaps = FALSE,
             quiet = FALSE)
### works only with matrices for the moment
{
    obj <- list(...)
    n <- length(obj)
    if (n == 1) return(obj[[1]])
    for (i in 1:n)
        if (!is.matrix(obj[[1]]))
            stop("the 'cbind' method for \"DNAbin\" accepts only matrices")
    NR <- unlist(lapply(obj, nrow))
    for (i in 1:n) class(obj[[i]]) <- NULL
    if (check.names) {
        nms <- unlist(lapply(obj, rownames))
        if (fill.with.gaps) {
            NC <- unlist(lapply(obj, ncol))
            nms <- unique(nms)
            ans <- matrix(as.raw(4), length(nms), sum(NC))
            rownames(ans) <- nms
            from <- 1
            for (i in 1:n) {
                to <- from + NC[i] - 1
                tmp <- rownames(obj[[i]])
                nmsi <- tmp[tmp %in% nms]
                ans[nmsi, from:to] <- obj[[i]][nmsi, , drop = FALSE]
                from <- to + 1
            }
        } else {
            tab <- table(nms)
            ubi <- tab == n
            nms <- names(tab)[which(ubi)]
            ans <- obj[[1]][nms, , drop = FALSE]
            for (i in 2:n)
                ans <- cbind(ans, obj[[i]][nms, , drop = FALSE])
            if (!quiet && !all(ubi))
                warning("some rows were dropped.")
        }
    } else {
        if (length(unique(NR)) > 1)
            stop("matrices do not have the same number of rows.")
        ans <- matrix(unlist(obj), NR)
        rownames(ans) <- rownames(obj[[1]])
    }
    class(ans) <- "DNAbin"
    ans
}

c.DNAbin <- function(..., recursive = FALSE)
{
    if (!all(unlist(lapply(list(...), is.list))))
        stop("the 'c' method for \"DNAbin\" accepts only lists")
    structure(NextMethod("c"), class = "DNAbin")
}

print.DNAbin <- function(x, printlen = 6, digits = 3, ...)
{
    if (is.list(x)) {
        n <- length(x)
        nms <- names(x)
        if (n == 1) {
            cat("1 DNA sequence in binary format stored in a list.\n\n")
            nTot <- length(x[[1]])
            cat("Sequence length:", nTot, "\n\n")
            cat("Label:", nms, "\n\n")
        } else {
            cat(n, "DNA sequences in binary format stored in a list.\n\n")
            tmp <- unlist(lapply(x, length))
            nTot <- sum(tmp)
            mini <- range(tmp)
            maxi <- mini[2]
            mini <- mini[1]
            if (mini == maxi)
                cat("All sequences of same length:", maxi, "\n")
            else {
                cat("Mean sequence length:", round(mean(tmp), 3), "\n")
                cat("   Shortest sequence:", mini, "\n")
                cat("    Longest sequence:", maxi, "\n")
            }
            TAIL <- "\n\n"
            if (printlen < n) {
                nms <- nms[1:printlen]
                TAIL <- "...\n\n"
            }
            cat("\nLabels:", paste(nms, collapse = " "), TAIL)
        }
    } else {
        nTot <- length(x)
        if (is.matrix(x)) {
            nd <- dim(x)
            nms <- rownames(x)
            cat(nd[1], "DNA sequences in binary format stored in a matrix.\n\n")
            cat("All sequences of same length:", nd[2], "\n")
            TAIL <- "\n\n"
            if (printlen < nd[1]) {
                nms <- nms[1:printlen]
                TAIL <- "...\n\n"
            }
            cat("\nLabels:", paste(nms, collapse = " "), TAIL)
        } else {
            cat("1 DNA sequence in binary format stored in a vector.\n\n")
            cat("Sequence length:", nTot, "\n\n")
        }
    }
    if (nTot <= 1e6) {
        cat("Base composition:\n")
        print(round(base.freq(x), digits))
    } else cat("More than 1 million nucleotides: not printing base composition\n")
}

as.DNAbin <- function(x, ...) UseMethod("as.DNAbin")

._cs_ <- c("a", "g", "c", "t", "r", "m", "w", "s", "k",
           "y", "v", "h", "d", "b", "n", "-", "?")

._bs_ <- c(136, 72, 40, 24, 192, 160, 144, 96, 80,
           48, 224, 176, 208, 112, 240, 4, 2)

## by Klaus:
as.DNAbin.character <- function(x, ...)
{
    ans <- as.raw(._bs_)[match(x, ._cs_)]
    if (is.matrix(x)) {
        dim(ans) <- dim(x)
        dimnames(ans) <- dimnames(x)
    }
    class(ans) <- "DNAbin"
    ans
}

as.DNAbin.alignment <- function(x, ...)
{
    n <- x$nb
    x$seq <- tolower(x$seq)
    ans <- matrix("", n, nchar(x$seq[1]))
    for (i in 1:n)
        ans[i, ] <- strsplit(x$seq[i], "")[[1]]
    rownames(ans) <- gsub(" +$", "", gsub("^ +", "", x$nam))
    as.DNAbin.character(ans)
}

as.DNAbin.list <- function(x, ...)
{
    obj <- lapply(x, as.DNAbin)
    class(obj) <- "DNAbin"
    obj
}

as.character.DNAbin <- function(x, ...)
{
    f <- function(xx) {
        ans <- ._cs_[match(as.numeric(xx), ._bs_)]
        if (is.matrix(xx)) {
            dim(ans) <- dim(xx)
            dimnames(ans) <- dimnames(xx)
        }
        ans
    }
    if (is.list(x)) lapply(x, f) else f(x)
}

base.freq <- function(x, freq = FALSE, all = FALSE)
{
    if (!inherits(x, "DNAbin")) stop('base.freq requires an object of class "DNAbin"')

    f <- function(x)
        .C(BaseProportion, x, as.integer(length(x)), double(17),
           NAOK = TRUE)[[3]]

    if (is.list(x)) {
        BF <- rowSums(sapply(x, f))
        n <- sum(sapply(x, length))
    } else {
        n <- length(x)
        BF <- .C(BaseProportion, x, as.integer(n), double(17),
                 NAOK = TRUE)[[3]]
    }

    names(BF) <- c("a", "c", "g", "t", "r", "m", "w", "s",
                   "k", "y", "v", "h", "d", "b", "n", "-", "?")
    if (all) {
        if (!freq) BF <- BF / n
    } else {
        BF <- BF[1:4]
        if (!freq) BF <- BF / sum(BF)
    }
    BF
}

Ftab <- function(x, y = NULL)
{
    if (is.null(y)) {
        if (is.list(x)) {
            y <- x[[2]]
            x <- x[[1]]
            if (length(x) != length(y))
                stop("'x' and 'y' not of the same length")
        } else { # 'x' is a matrix
            y <- x[2, , drop = TRUE]
            x <- x[1, , drop = TRUE]
        }
    } else {
        x <- as.vector(x)
        y <- as.vector(y)
        if (length(x) != length(y))
            stop("'x' and 'y' not of the same length")
    }
    out <- matrix(0, 4, 4)
    k <- c(136, 40, 72, 24)
    for (i in 1:4) {
        a <- x == k[i]
        for (j in 1:4) {
            b <- y == k[j]
            out[i, j] <- sum(a & b)
        }
    }
    dimnames(out)[1:2] <- list(c("a", "c", "g", "t"))
    out
}

GC.content <- function(x) sum(base.freq(x)[2:3])

seg.sites <- function(x)
{
    if (is.list(x)) x <- as.matrix(x)
    if (is.vector(x)) n <- 1
    else { # 'x' is a matrix
        n <- dim(x)
        s <- n[2]
        n <- n[1]
    }
    if (n == 1) return(integer(0))
    ans <- .C(SegSites, x, as.integer(n), as.integer(s),
              integer(s), NAOK = TRUE)
    which(as.logical(ans[[4]]))
}

dist.dna <- function(x, model = "K80", variance = FALSE, gamma = FALSE,
                     pairwise.deletion = FALSE, base.freq = NULL,
                     as.matrix = FALSE)
{
    MODELS <- c("RAW", "JC69", "K80", "F81", "K81", "F84", "T92", "TN93",
                "GG95", "LOGDET", "BH87", "PARALIN", "N", "TS", "TV",
                "INDEL", "INDELBLOCK")
    imod <- pmatch(toupper(model), MODELS)
    if (is.na(imod))
        stop(paste("'model' must be one of:",
                   paste("\"", MODELS, "\"", sep = "", collapse = " ")))
    if (imod == 11 && variance) {
        warning("computing variance not available for model BH87")
        variance <- FALSE
    }
    if (gamma && imod %in% c(1, 5:7, 9:17)) {
        warning(paste("gamma-correction not available for model", model))
        gamma <- FALSE
    }
    if (is.list(x)) x <- as.matrix(x)
    nms <- dimnames(x)[[1]]
    n <- dim(x)
    s <- n[2]
    n <- n[1]

    if (imod %in% c(4, 6:8)) {
        BF <- if (is.null(base.freq)) base.freq(x) else base.freq
    } else BF <- 0

    if (imod %in% 16:17) pairwise.deletion <- TRUE

    if (!pairwise.deletion) {
        keep <- .C(GlobalDeletionDNA, x, n, s, rep(1L, s))[[4]]
        x <- x[, as.logical(keep)]
        s <- dim(x)[2]
    }
    Ndist <- if (imod == 11) n*n else n*(n - 1)/2
    var <- if (variance) double(Ndist) else 0
    if (!gamma) gamma <- alpha <- 0
    else {
        alpha <- gamma
        gamma <- 1
    }
    d <- .C(dist_dna, x, as.integer(n), as.integer(s), imod,
            double(Ndist), BF, as.integer(pairwise.deletion),
            as.integer(variance), var, as.integer(gamma),
            as.double(alpha), NAOK = TRUE)
    if (variance) var <- d[[9]]
    d <- d[[5]]
    if (imod == 11) {
        dim(d) <- c(n, n)
        dimnames(d) <- list(nms, nms)
    } else {
        attr(d, "Size") <- n
        attr(d, "Labels") <- nms
        attr(d, "Diag") <- attr(d, "Upper") <- FALSE
        attr(d, "call") <- match.call()
        attr(d, "method") <- model
        class(d) <- "dist"
        if (as.matrix) d <- as.matrix(d)
    }
    if (variance) attr(d, "variance") <- var
    d
}

image.DNAbin <- function(x, what, col, bg = "white", xlab = "", ylab = "",
                         show.labels = TRUE, cex.lab = 1, legend = TRUE, ...)
{
    what <-
        if (missing(what)) c("a", "g", "c", "t", "n", "-") else tolower(what)
    if (missing(col))
        col <- c("red", "yellow", "green", "blue", "grey", "black")
    n <- (dx <- dim(x))[1] # number of sequences
    s <- dx[2] # number of sites
    y <- integer(N <- length(x))
    ncl <- length(what)
    col <- rep(col, length.out = ncl)
    brks <- 0.5:(ncl + 0.5)
    sm <- 0L
    for (i in ncl:1) {
        k <- ._bs_[._cs_ == what[i]]
        sel <- which(x == k)
        if (L <- length(sel)) {
            y[sel] <- i
            sm <- sm + L
        } else {
            what <- what[-i]
            col <- col[-i]
            brks <- brks[-i]
        }
    }
    dim(y) <- dx
    ## if there's no 0 in y, must drop 'bg' from the cols passed to image:
    if (sm == N) {
        leg.co <- co <- col
        leg.txt <- toupper(what)
    } else {
        co <- c(bg, col)
        leg.txt <- c(toupper(what), "others")
        leg.co <- c(col, bg)
        brks <- c(-0.5, brks)
    }
    yaxt <- if (show.labels) "n" else "s"
    image.default(1:s, 1:n, t(y), col = co, xlab = xlab,
                  ylab = ylab, yaxt = yaxt, breaks = brks, ...)
    if (show.labels)
        mtext(rownames(x), side = 2, line = 0.1, at = 1:n,
              cex = cex.lab, adj = 1, las = 1)
    if (legend) {
        psr <- par("usr")
        xx <- psr[2]/2
        yy <- psr[4] * (0.5 + 0.5/par("plt")[4])
        legend(xx, yy, legend = leg.txt, pch = 22, pt.bg = leg.co,
               pt.cex = 2, bty = "n", xjust = 0.5, yjust = 0.5,
               horiz = TRUE, xpd = TRUE)
    }
}

alview <- function(x, file = "", uppercase = TRUE)
{
    if (is.list(x)) x <- as.matrix(x)
    taxa <- formatC(labels(x), width = -1)
    x <- as.character(x)
    s <- ncol(x)
    if (nrow(x) > 1) {
        for (j in seq_len(s)) {
            q <- which(x[-1L, j] == x[1L, j]) + 1L
            x[q, j] <- "."
        }
    }
    x <- apply(x, 1L, paste, collapse = "")
    if (uppercase) x <- toupper(x)
    res <- paste(taxa, x)
    hdr <- formatC(s, width = nchar(res[1]))
    cat(hdr, file = file, sep = "\n")
    cat(res, file = file, sep = "\n", append = TRUE)
}

where <- function(x, pattern)
{
    pat <- as.DNAbin(strsplit(pattern, NULL)[[1]])
    p <- as.integer(length(pat))

    foo <- function(x, pat, p) {
        s <- as.integer(length(x))
        if (s < p) stop("sequence shorter than the pattern")
        ans <- .C(C_where, x, pat, s, p, integer(s), 0L, NAOK = TRUE)
        n <- ans[[6]]
        if (n) ans[[5]][seq_len(n)] - p + 2L else integer()
    }

    if (is.list(x)) return(lapply(x, foo, pat = pat, p = p))
    if (is.matrix(x)) {
        n <- nrow(x)
        res <- vector("list", n)
        for (i in seq_len(n))
            res[[i]] <- foo(x[i, , drop = TRUE], pat, p)
        names(res) <- rownames(x)
        return(res)
    }
    foo(x, pat, p) # if x is a vector
}

## conversions from BioConductor:

.DNAString2DNAbin <- function(from)
    as.DNAbin.character(strsplit(tolower(as.character(from)), "")[[1]])

.DNAStringSet2DNAbin <- function(from)
    structure(lapply(from, .DNAString2DNAbin), class = "DNAbin")

.DNAAlignment2DNAbinMatrix <- function(from)
    as.matrix(.DNAStringSet2DNAbin(as(from, "DNAStringSet")))

## S3 versions:
as.DNAbin.DNAString <- function(x, ...) .DNAString2DNAbin(x)
as.DNAbin.DNAStringSet <- function(x, ...) .DNAStringSet2DNAbin(x)
as.DNAbin.PairwiseAlignmentsSingleSubject <- function(x, ...) .DNAAlignment2DNAbinMatrix(x)
as.DNAbin.DNAMultipleAlignment <- function(x, ...) .DNAAlignment2DNAbinMatrix(x)

## S4 versions:
## setAs("DNAString", "DNAbin", .DNAString2DNAbin)
## setAs("DNAStringSet", "DNAbin", .DNAStringSet2DNAbin)
## setAs("PairwiseAlignmentsSingleSubject", "DNAbin", .DNAAlignment2DNAbinMatrix)
## setAs("DNAMultipleAlignment", "DNAbin", .DNAAlignment2DNAbinMatrix)
