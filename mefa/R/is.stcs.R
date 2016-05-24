`is.stcs` <-
function(x) {
    inherits(x, "stcs") && colnames(x) == c("samp", "taxa", "count", "segm") && is.numeric(x[, 3]) &&
    is.factor(x$samp) && is.factor(x$taxa) && is.factor(x$segm)
}

