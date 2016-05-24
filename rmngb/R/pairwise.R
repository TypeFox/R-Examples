pairwise.chisq.test <- function(x, ...)
    UseMethod("pairwise.chisq.test")

pairwise.chisq.test.table <- function(x, p.adj = p.adjust.methods,
                                      DNAME = NULL, ...) {
    if (is.null(DNAME))
        DNAME <- deparse(substitute(tab))
    p.adj <- match.arg(p.adj)
    
    k <- dim(x)[1]
    p.value <- rep(NA, k ^ 2)
    for (row1 in 1:(k - 1)) {
        for (row2 in (row1 + 1):k) {
            xi <- asInteger(k * (row1 - 1) + row2)
            p.value[xi] <- chisq.test(x[c(row1, row2), ])$p.value
        }
    }
    
    p.value <- p.adjust(p.value, method = p.adj)
    dn <- list(dimnames(x)[[1]],
               dimnames(x)[[1]])
    p.value <- matrix(p.value, nrow = k, dimnames = dn)
    
    structure(list(method = "Pearson's Chi-squared tests",
                   data.name = DNAME,
                   p.value = p.value,
                   p.adjust.method = p.adj),
              class = "pairwise.htest")
}

pairwise.chisq.test.default <- function(x, g, p.adj = p.adjust.methods, ...) {
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    p.adj <- match.arg(p.adj)
    tab <- table(g, x)
    
    pairwise.chisq.test(tab, p.adj = p.adj, DNAME = DNAME)
}

pairwise.fisher.test <- function(x, ...)
    UseMethod("pairwise.fisher.test")

pairwise.fisher.test.table <- function(x, p.adj = p.adjust.methods,
                                       DNAME = NULL, ...) {
    if (is.null(DNAME))
        DNAME <- deparse(substitute(tab))
    p.adj <- match.arg(p.adj)
    
    k <- dim(x)[1]
    p.value <- rep(NA, k ^ 2)
    for (row1 in 1:(k - 1)) {
        for (row2 in (row1 + 1):k) {
            xi <- asInteger(k * (row1 - 1) + row2)
            p.value[xi] <- fisher.test(x[c(row1, row2), ])$p.value
        }
    }
    
    p.value <- p.adjust(p.value, method = p.adj)
    dn <- list(dimnames(x)[[1]],
               dimnames(x)[[1]])
    p.value <- matrix(p.value, nrow = k, dimnames = dn)
    
    structure(list(method = "Fisher's Exact Tests for Count Data",
                   data.name = DNAME,
                   p.value = p.value,
                   p.adjust.method = p.adj),
              class = "pairwise.htest")
}

pairwise.fisher.test.default <- function(x, g, p.adj = p.adjust.methods, ...) {
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    p.adj <- match.arg(p.adj)
    tab <- table(g, x)
    
    pairwise.fisher.test(tab, p.adj = p.adj, DNAME = DNAME)
}
