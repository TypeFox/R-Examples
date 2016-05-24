############################################################################################
## package 'secrlinear'
## sim.linearpopn.R
## last changed 2014-08-29, 2014-09-03
############################################################################################

sim.linearpopn <- 	function (mask, D, N, Ndist = c('poisson', 'fixed'), ...) {
    Ndist <- match.arg(Ndist)
    if (!inherits(mask, 'linearmask'))
        stop ("requires linearmask")
    if (missing(D)) {
        D <- rep(N / masklength(mask), nrow(mask))
    }
    sim.popn(core = mask, D = D, model2D = "linear", Ndist = Ndist, ...)
}
