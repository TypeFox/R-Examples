## S4 -> S3

as.mefa.Mefa <- function(x, ...) {
    mefa(xtab = x@xtab, samp = x@samp, taxa = x@taxa, 
        id.samp = NULL, id.taxa = NULL,
        segment = FALSE, nested = FALSE, drop.zero = FALSE, drop.index = FALSE,
        xtab.fixed = x@join == "left")
}

as.mefa.sparseMatrix <- function(x, ...) {
    mefa(xtab = x, samp = NULL, taxa = NULL, 
        id.samp = NULL, id.taxa = NULL,
        segment = FALSE, nested = FALSE, drop.zero = FALSE, drop.index = FALSE,
        xtab.fixed = TRUE)
}

as.stcs.sparseMatrix <- function(x, ...) {
    melt(as.mefa(x))
}

as.stcs.Mefa <- function(x, ...) {
    melt(as.mefa(x))
}

melt.sparseMatrix <- as.stcs.sparseMatrix

melt.Mefa <- as.stcs.Mefa

## S3 -> S4

as.Xtab <- function (x, ...) 
    UseMethod("as.Xtab")

as.Mefa <- function (x, ...) 
    UseMethod("as.Mefa")

as.Xtab.mefa <- function(x, ...) {
    as(x$xtab, "dgCMatrix")
}

as.Mefa.mefa <- function(x, ...) {
    new("Mefa", 
        xtab = as(x$xtab, "dgCMatrix"), 
        samp = x$samp,
        taxa = x$taxa,
        join = ifelse(attr(x, "xtab.fixed"), "left", "inner"))
}

as.Xtab.stcs <- function(x, ...) {
    mefa4::Xtab(count ~ samp + taxa, x)
}

as.Mefa.stcs <- function(x, ...) {
    mefa4::Mefa(as.Xtab(x))
}
