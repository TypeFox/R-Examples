GHpoints <-
function (form, k) {
    av <- all.vars(form)
    factors <- sum(av %in% c("z1", "z2"))
    cf <- paste(form[3])
    #form. <- paste(" ~ ", paste("z", 1:factors, collapse = " + ", sep = ""), " + ", cf)
    form. <- paste(" ~ ", cf)
    form. <- as.formula(form.)
    GH <- gauher(k)
    grid.t <- expand.grid(lapply(1:factors, function (k, u) u$x, u = GH))
    names(grid.t) <- paste("z", 1:factors, sep = "")
    out <- model.matrix(form., sqrt(2) * grid.t)
    colnams <- colnames(out)
    dimnames(out) <- attr(out, "assign") <- NULL
    grid.w <- as.matrix(expand.grid(lapply(1:factors, function (k, u) u$w, u = GH)))
    grid.w <- (2^(factors/2)) * apply(grid.w, 1, prod) * exp(rowSums(grid.t * grid.t))
    grid.w <- grid.w * exp(rowSums(dnorm(out[, seq(2, factors + 1), drop = FALSE], log = TRUE)))
    names(grid.w) <- NULL
    list(x = out, w = grid.w, colnams = colnams)
}
