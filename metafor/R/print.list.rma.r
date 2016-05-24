print.list.rma <-
function (x, digits, ...) 
{
    if (!is.element("list.rma", class(x))) 
        stop("Argument 'x' must be an object of class \"list.rma\".")
    if (missing(digits)) 
        digits <- x$digits
    attr(x, "class") <- NULL
    slab.pos <- which(names(x) == "slab")
    out <- x[seq_len(slab.pos - 1)]
    out <- data.frame(out, row.names = x$slab)
    if (nrow(out) == 0L) 
        stop("All values are NA.", call. = FALSE)
    transf.true <- 0
    if (exists("transf", where = x, inherits = FALSE) && x$transf) {
        transf.true <- 1
        out$se <- NULL
    }
    if (exists("method", where = x, inherits = FALSE)) {
        min.pos <- slab.pos - is.element("tau2.level", names(x)) - 
            is.element("gamma2.level", names(x)) - is.element("X", 
            names(x)) - transf.true
    }
    else {
        min.pos <- slab.pos - transf.true
    }
    sav <- out[, seq_len(min.pos - 1)]
    out[, seq_len(min.pos - 1)] <- apply(out[, seq_len(min.pos - 
        1), drop = FALSE], 2, formatC, digits = digits, format = "f")
    print(out, quote = FALSE, right = TRUE)
    invisible(sav)
}
