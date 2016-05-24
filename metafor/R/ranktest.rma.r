ranktest.rma <-
function (x, ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    if (is.element("robust.rma", class(x))) 
        stop("Function not applicable to objects of class \"robust.rma\".")
    yi <- x$yi
    vi <- x$vi
    res <- rma(yi, vi, method = "FE")
    b <- res$b
    vb <- res$vb
    vi.star <- vi - c(vb)
    yi.star <- (yi - c(b))/sqrt(vi.star)
    res <- cor.test(yi.star, vi, method = "kendall", exact = TRUE)
    pval <- res$p.value
    tau <- c(res$estimate)
    res <- list(tau = tau, pval = pval, digits = x$digits)
    class(res) <- c("ranktest.rma")
    return(res)
}
