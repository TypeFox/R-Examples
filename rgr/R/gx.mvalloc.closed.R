gx.mvalloc.closed <-
function (pcrit = 0.05, xx, ...) 
{
    if (!is.matrix(xx)) 
        stop("  ", deparse(substitute(xx)), " is not a Matrix")
    temp.x <- remove.na(xx)
    x <- clr(temp.x$x)
    nx <- temp.x$n
    px <- temp.x$m
    ListOfGroups <- list(...)
    kk <- length(ListOfGroups)
    cat("  k =", kk, "\t px =", px, "\t nx =", nx, "\n")
    temp <- (nx - px)/(px * (nx + 1))
    groups <- character(kk)
    md <- numeric(nx * kk)
    md <- array(md, dim = c(nx, kk))
    pgm <- numeric(nx * kk)
    pgm <- array(pgm, dim = c(nx, kk))
    xalloc <- integer(nx)
    for (k in 1:kk) {
        if (ListOfGroups[[k]]$p != px) 
            stop("\n  p != px for data set ", k, "\n")
        groups[k] <- ListOfGroups[[k]]$main
        if (ListOfGroups[[k]]$nc <= 5 * ListOfGroups[[k]]$p) 
            cat("  *** Proceed with Care, n < 5p for group", 
                groups[k], "***\n")
        if (ListOfGroups[[k]]$nc <= 3 * ListOfGroups[[k]]$p) 
            cat("  *** Proceed with Great Care, n < 3p for group", 
                groups[k], "***\n")
        md[, k] <- mahalanobis(x, ListOfGroups[[k]]$mean, ListOfGroups[[k]]$cov.inv, 
            inverted = TRUE)
        pgm[, k] <- round(1 - pf(temp * md[, k], px, nx - px), 
            4)
        md[, k] <- md[, k] + det(ListOfGroups[[k]]$cov)
    }
    for (i in 1:nx) {
        xalloc[i] <- which(md[i, ] == min(md[i, ]))
        if (max(pgm[i, ]) < pcrit) 
            xalloc[i] <- 0
    }
    invisible(list(groups = groups, kk = kk, n = nx, p = px, 
        pcrit = pcrit, pgm = pgm, xalloc = xalloc))
}
