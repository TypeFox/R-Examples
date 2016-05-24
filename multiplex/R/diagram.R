diagram <-
function (x, unord = TRUE, attrs = NULL, main = NULL, cex.main = graphics::par()$cex.main) 
{
    if (requireNamespace("Rgraphviz", quietly = TRUE)) {
        if (is.null(dimnames(x)[[1]]) == TRUE) 
            rownames(x) <- colnames(x) <- as.character(utils::as.roman(c(1:dim(x)[1])))
        if (is.null(attrs) == TRUE) 
            attrs = list(graph = list(rankdir = "BT"), edge = list(arrowsize = "0", 
                minlen = "1"), node = list(shape = "rectangle", 
                color = "white", fixedsize = FALSE))
        po <- x & (1L - t(x))
        diag(po) <- 0L
        if (isTRUE(ncol(x) == nrow(x)) == TRUE) {
            for (i in seq_len(ncol(po))) {
                tmp <- outer(po[, i], po[i, ], pmin.int)
                po <- pmin(po, (1L - tmp))
            }
            rm(tmp)
        }
        else {
            stop("binary operation on non-conformable arrays")
        }
        if (unord == FALSE) {
            px <- po
            out <- vector()
            k <- 1L
            for (i in 1:nrow(px)) {
                if (sum(px[i, ] + px[, i]) == 0L) {
                  out[k] <- i
                  k <- k + 1L
                }
            }
            rm(i)
            d <- nrow(px)
            for (j in out) {
                for (k in 1:d) {
                  px[, j] <- px[j, ] <- NA
                }
            }
            rm(j)
            lb <- dimnames(px)[[1]]
            for (l in out) {
                lb[l] <- NA
            }
            rm(l)
            npx <- data.frame(matrix(0L, ncol = (nrow(px) - length(out)), 
                nrow = 0L))
            colnames(npx) <- as.vector(stats::na.exclude(lb))
            for (i in 1:d) {
                ifelse(isTRUE(all(is.na(px[i, ])) == FALSE) == 
                  TRUE, npx[i, ] <- as.vector(stats::na.exclude(px[i, 
                  ])), NA)
            }
            rm(i)
            po <- as.matrix(stats::na.exclude(npx))
            dimnames(po)[[1]] <- dimnames(po)[[2]] <- as.vector(stats::na.exclude(lb))
        }
        Rgraphviz::plot(methods::as(po, "graphNEL"), attrs = attrs, 
            main = main, cex.main = cex.main)
    }
    else stop("Package 'Rgraphviz' needs to be properly installed")
}
