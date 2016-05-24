varianceDecomposition <-
function (obj) 
{
    ans <- list()
    if (class(obj) == "noia.linear" || class(obj) == "noia.linear.gpmap" || 
        class(obj) == "noia.multilinear") {
        n <- names(obj$variances)
        if (obj$nloc == 1) {
            n.a <- c(0, 1, 0)
            n.d <- c(0, 0, 1)
            n.e <- c(0, 0, 0)
        }
        else {
            n.a <- apply(sapply(strsplit(n, ""), "c") == noia::effectsNames[2], 
                2, "sum")
            n.d <- apply(sapply(strsplit(n, ""), "c") == noia::effectsNames[3], 
                2, "sum")
            n.e <- apply(sapply(strsplit(n, ""), "c") == noia::effectsNames[4], 
                2, "sum")
        }
        for (lev in 1:(obj$nloc)) {
            if (class(obj) == "noia.linear" || class(obj) == 
                "noia.linear.gpmap" || lev < 2) {
                for (nr.d in 0:lev) {
                  nr.a <- lev - nr.d
                  v <- sum(obj$variances[(n.a == nr.a) & (n.d == 
                    nr.d)])
                  {
                    order.label <- as.character(lev)
                    component.label <- paste(c(rep("A", nr.a), 
                      rep("D", nr.d)), collapse = "")
                    names(v) <- component.label
                    ans[[order.label]] <- c(ans[[order.label]], 
                      v)
                  }
                }
            }
            else {
                v <- sum(obj$variances[n.e == lev])
                if (is.finite(v)) {
                  order.label <- as.character(lev)
                  component.label <- paste(rep("E", lev), collapse = "")
                  names(v) <- component.label
                  ans[[order.label]] <- c(ans[[order.label]], 
                    v)
                }
            }
        }
    }
    else {
        stop("Class", class(obj), "unknown.\n")
    }
    if (class(obj) == "noia.linear.gpmap") {
        ans$V_G <- obj$V_G
    }
    class(ans) <- c("noia.vardec", class(ans))
    return(ans)
}
