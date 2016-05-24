HPDinterval <- function(obj, prob = 0.95, ...) UseMethod("HPDinterval")

HPDinterval.mcmc <- function(obj, prob = 0.95, ...)
{
    obj <- as.matrix(obj)
    vals <- apply(obj, 2, sort)
    if (!is.matrix(vals)) stop("obj must have nsamp > 1")
    nsamp <- nrow(vals)
    npar <- ncol(vals)
    gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
    init <- 1:(nsamp - gap)
    inds <- apply(vals[init + gap, ,drop=FALSE] - vals[init, ,drop=FALSE],
                  2, which.min)
    ans <- cbind(vals[cbind(inds, 1:npar)],
                 vals[cbind(inds + gap, 1:npar)])
    dimnames(ans) <- list(colnames(obj), c("lower", "upper"))
    attr(ans, "Probability") <- gap/nsamp
    ans
}

HPDinterval.mcmc.list <- function(obj, prob = 0.95, ...)
    lapply(obj, HPDinterval, prob)
