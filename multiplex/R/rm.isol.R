rm.isol <-
function (x, diag.incl = TRUE) 
{
    ifelse(isTRUE(all(x == 0)) == TRUE, return(NULL), NA)
    if (isTRUE(dim(x)[3] > 1L) == TRUE) {
        tmp <- x
        ifelse(isTRUE(is.null(dimnames(x)[[1]])) == TRUE & isTRUE(is.null(dimnames(x)[[2]])) == 
            TRUE, dimnames(x)[[2]] <- dimnames(x)[[1]] <- 1:dim(x)[1], 
            NA)
        mpx <- data.frame(matrix(0L, ncol = ncol(x), nrow = nrow(x)))
        mpx <- x[, , 1] + x[, , 2]
        if (isTRUE(dim(x)[3] > 2L) == TRUE) {
            for (m in 3:dim(x)[3]) {
                mpx <- mpx + x[, , m]
            }
            rm(m)
        }
        if (isTRUE(diag.incl == FALSE) == TRUE) 
            diag(mpx) <- 0L
        px <- mpx
        rm(mpx)
        out <- vector()
        for (i in 1:nrow(px)) {
            if (isTRUE(sum(px[i, ] + px[, i]) == 0L) == TRUE) {
                out[length(out) + 1] <- i
            }
        }
        rm(i)
        for (j in out) {
            for (k in 1:nrow(px)) {
                px[, j] <- px[j, ] <- NA
            }
        }
        rm(j)
        if (is.null(dimnames(x)) == TRUE) 
            dimnames(px)[[1]] <- 1:dim(x)[1]
        lb <- dimnames(px)[[1]]
        for (l in out) {
            lb[l] <- NA
        }
        rm(l)
        pxx <- x
        for (m in 1:dim(x)[3]) {
            pxx[which(is.na(lb)), , m] <- NA
            pxx[, which(is.na(lb)), m] <- NA
        }
        rm(m)
        mx <- array(dim = c((nrow(px) - length(out)), (nrow(px) - 
            length(out)), dim(x)[3]))
        for (m in 1:dim(x)[3]) {
            npx <- data.frame(matrix(0L, ncol = (nrow(px) - length(out)), 
                nrow = 0L))
            colnames(npx) <- as.vector(stats::na.exclude(lb))
            for (i in 1:nrow(px)) {
                ifelse(isTRUE(all(is.na(pxx[i, , m])) == FALSE) == 
                  TRUE, npx[i, ] <- as.vector(stats::na.exclude(pxx[i, 
                  , m])), NA)
            }
            rm(i)
            mx[, , m] <- as.matrix(stats::na.exclude(npx))
        }
        rm(m)
        if (is.null(dimnames(x)) == FALSE) 
            dimnames(mx)[[1]] <- dimnames(mx)[[2]] <- as.vector(stats::na.exclude(lb))
        if (isTRUE(is.null(dimnames(tmp)[[1]])) == TRUE & isTRUE(is.null(dimnames(tmp)[[2]])) == 
            TRUE) {
            dimnames(mx) <- NULL
            warning("No labels are provided in 'x'")
        }
        if (is.null(dimnames(x)[[3]]) == FALSE) 
            dimnames(mx)[[3]] <- dimnames(x)[[3]]
        return(mx)
    }
    else {
        if (isTRUE(diag.incl == FALSE) == TRUE) {
            dg <- diag
            diag(x) <- 0L
        }
        if (isTRUE(sum(x) == 0L) == TRUE) {
            x
        }
        else {
            px <- x
            out <- vector()
            for (i in 1:nrow(px)) {
                if (isTRUE(sum(px[i, ] + px[, i]) == 0L) == TRUE) {
                  out[length(out) + 1] <- i
                }
            }
            rm(i)
            for (j in out) {
                for (k in 1:nrow(px)) {
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
            for (i in 1:nrow(px)) {
                ifelse(isTRUE(all(is.na(px[i, ])) == FALSE) == 
                  TRUE, npx[i, ] <- as.vector(stats::na.exclude(px[i, 
                  ])), NA)
            }
            rm(i)
            mx <- as.matrix(stats::na.exclude(npx))
            dimnames(mx)[[1]] <- dimnames(mx)[[2]] <- as.vector(stats::na.exclude(lb))
            mx
        }
    }
}
