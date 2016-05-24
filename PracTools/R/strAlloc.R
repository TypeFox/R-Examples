strAlloc <- function (n.tot = NULL, Nh = NULL, Sh = NULL, cost = NULL, ch = NULL,
    V0 = NULL, CV0 = NULL, ybarU = NULL, alloc)
{
    if (!(alloc %in% c("prop", "neyman", "totcost", "totvar")))
        stop("Illegal allocation specified.\n")
    if (is.null(Nh))
        stop("Nh cannot be NULL.\n")
    if (alloc == "prop") {
        if (is.null(n.tot))
            stop("n.tot is required for proportional allocation.\n")
        if (!is.null(Sh) || !is.null(cost) || !is.null(ch))
            warning("Sh, cost, and ch are ignored for proportional allocation.\n")
    }
    if (!(alloc == "prop")) {
        if (is.null(Sh))
            stop("Sh cannot be NULL unless allocation is proportional.\n")
        if (!(alloc == "neyman") & !(alloc == "totcost")) {
            if (!is.null(CV0) & !is.null(V0))
                stop("Only one of CV0 and V0 should be non-null unless allocation is proportional, Neyman, cost 
                constrained.\n")
        }
    }
    if (alloc == "neyman") {
        if (is.null(n.tot)) {
            stop("n.tot is required for Neyman allocation.\n")
        }
        if (!is.null(cost)) {
            warning("Cost is ignored when allocation is Neyman.\n")
        }
    }
    if ((alloc %in% c("totvar", "totcost")) & !is.null(n.tot)) {
        stop("n.tot must be NULL if allocation is totvar or totcost.\n")
    }
    if (any(Nh <= 0, Sh <= 0, CV0 <= 0, V0 <= 0, cost <= 0))
        stop("Nh, Sh, CV0, V0, and cost cannot be <= 0.\n")
    if (alloc == "totcost") {
        if (sum(ch) > cost)
            stop("Sum of stratum unit costs cannot exceed total cost.\n")
    }
    N <- sum(Nh)
    Wh <- Nh/N
    if (alloc == "prop")
        nh <- n.tot * Wh
    if ((alloc == "neyman") & !is.null(n.tot)) {
        nh <- n.tot * Wh * Sh/sum(Wh * Sh)
    }
    if ((alloc == "totcost")) {
        if (is.null(cost))
            stop("If alloc=totcost, cost must be specified.\n")
        d1 <- sum(Wh * Sh/sqrt(ch))
        ph.cost <- Wh * Sh/sqrt(ch)/d1
        n.cost <- cost * d1/sum(Wh * Sh * sqrt(ch))
        nh <- n.cost * ph.cost
    }
    if ((alloc == "totvar")) {
        if (is.null(CV0) & is.null(V0))
            stop("CV0 and V0 cannot both be NULL if allocation is totvar.\n")
        if (!is.null(CV0)) {
            V0 <- (CV0 * ybarU)^2
        }
        d1 <- sum(Wh * Sh/sqrt(ch))
        d2 <- sum(Wh * Sh * sqrt(ch))
        d3 <- V0 + sum(Wh * Sh^2)/N
        ph.cost <- Wh * Sh/sqrt(ch)/d1
        n.cost <- d1 * d2/d3
        nh <- n.cost * ph.cost
    }
    if (any(nh > Nh)) {
        viol <- (1:length(nh))[nh > Nh]
        warning("nh is larger than Nh in strata: ", viol, "\n")
    }
    if (alloc != "prop") {
        a.se <- sqrt(sum(Nh * (Nh/nh - 1) * (Sh^2))/N^2)
        structure(list(allocation = alloc, Nh = Nh, Sh = Sh, nh = nh, `nh/n` = nh/sum(nh),
            `anticipated SE of estimated mean` = a.se),
            class = "power.htest")
    }
    else {
        structure(list(allocation = alloc, Nh = Nh, nh = nh, `nh/n` = nh/sum(nh)), class = "power.htest")
    }
}
