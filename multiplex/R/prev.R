prev <-
function (x) 
{
    if (is.array(x) == FALSE) 
        stop("Data must be a stacked array of square matrices.")
    if (is.na(dim(x)[3]) == TRUE) {
        s0 <- data.frame(matrix(ncol = 1L, nrow = 1L))
        if (isTRUE(all.equal(replace(x %*% x, x %*% x >= 1L, 
            1L), x) == TRUE)) 
            s0[1, 1] <- 1L
        Bx <- array(dim = c(dim(x)[1], dim(x)[2], 2L))
        Bx[, , 1] <- as.matrix(x)
        Bx[, , 2] <- replace(x %*% x, x %*% x >= 1L, 1L)
    }
    if (is.na(dim(x)[3]) == FALSE) {
        tmp0 <- data.frame(matrix(ncol = (dim(x)[1] * dim(x)[2]), 
            nrow = 0L))
        for (i in 1:dim(x)[3]) {
            ifelse(isTRUE(dim(x)[3] > 1L) == TRUE, tmp0[i, ] <- as.vector(x[, 
                , i]), tmp0 <- as.vector(x))
        }
        rm(i)
        if (isTRUE(is.null(dim(tmp0)) == FALSE) == TRUE) 
            rownames(tmp0) <- dimnames(x)[[3]]
        if (isTRUE(dim(x)[3] < 2L) == TRUE) 
            x <- array(tmp0, c(dim(x)[1], dim(x)[2]))
        if (isTRUE(dim(x)[3] > 1L) == TRUE) {
            tmp <- array(dim = c(dim(x)[1], dim(x)[2], nrow(unique(tmp0))))
            for (i in 1:nrow(unique(tmp0))) {
                tmp[, , i][1:(dim(x)[1] * dim(x)[2])] <- as.numeric(unique(tmp0)[i, 
                  ])
            }
            rm(i)
            if (is.null(dimnames(tmp)[[1]]) == FALSE) 
                dimnames(tmp)[[3]] <- rownames(unique(tmp0))
            if (is.null(dimnames(x)[[1]]) == FALSE) 
                dimnames(tmp)[[1]] <- dimnames(tmp)[[2]] <- dimnames(x)[[1]]
            x <- tmp
            dimnames(x)[[3]] <- as.list(rownames(unique(tmp0)))
        }
        rm(tmp0, tmp)
        s0 <- data.frame(matrix(ncol = dim(x)[3], nrow = dim(x)[3]))
        for (k in 1:dim(x)[3]) {
            for (j in 1:dim(x)[3]) {
                tmp <- x[, , j] %*% x[, , k]
                tmp <- replace(tmp, tmp >= 1L, 1L)
                for (i in dim(x)[3]:1) {
                  if (isTRUE(all.equal(tmp, x[, , i]) == TRUE)) 
                    s0[j, k] <- i
                }
            }
        }
        rm(i, j, k)
        dimnames(s0)[[1]] <- 1:dim(x)[3]
        dimnames(s0)[[2]] <- 1:dim(x)[3]
        if (sum(as.numeric(is.na(s0))) == 0L) 
            Bx <- x
        if (sum(as.numeric(is.na(s0))) > 0L) {
            Bx <- array(dim = c(dim(x)[1], dim(x)[2], 0L))
            for (i in 1:nrow(s0)) {
                for (j in 1:length(which(is.na(s0[i, ])))) {
                  if (length(which(is.na(s0[i, ]))) > 0L) 
                    Bx <- zbnd(Bx, (replace(x[, , i] %*% x[, 
                      , which(is.na(s0[i, ]))[j]], x[, , i] %*% 
                      x[, , which(is.na(s0[i, ]))[j]] >= 1L, 
                      1L)))
                }
            }
            rm(i, j)
            tmp <- data.frame(matrix(ncol = (dim(x)[1] * dim(x)[2]), 
                nrow = 0L))
            for (i in 1:dim(Bx)[3]) {
                tmp[i, ] <- as.vector(Bx[, , i])
            }
            rm(i)
            xBx <- array(dim = c(dim(x)[1], dim(x)[2], nrow(unique(tmp))))
            for (i in 1:nrow(unique(tmp))) {
                xBx[, , i][1:(dim(Bx)[1] * dim(Bx)[2])] <- as.numeric(unique(tmp)[i, 
                  ])
            }
            rm(i)
            if (is.null(dimnames(xBx)) == FALSE) 
                dimnames(xBx)[[3]] <- (dim(x)[3] + 1L):(dim(xBx)[3] + 
                  dim(x)[3])
            Bx <- zbnd(x, xBx)
            rm(xBx, tmp)
        }
    }
    if (is.null(dimnames(x)[[3]]) == FALSE) 
        dimnames(s0)[[2]] <- dimnames(x)[[3]]
    pct <- round(length(attr(stats::na.omit(as.vector(unlist(s0))), 
        "na.action"))/length(as.vector(unlist(s0))), 2)
    d <- as.numeric(sort(unlist(s0), decreasing = TRUE))[1]
    if (isTRUE(d > 7L) == TRUE) {
        if (isTRUE(pct < 0.5) == TRUE) 
            return(list(`2stpT` = s0, PcU2stpT = pct, ordr = d))
        if (isTRUE(pct > 0.5) == TRUE) 
            return(list(`2stpT` = s0, PcU2stpT = pct, ordr = d, 
                Note = c("Complete semigroup construction may take long time")))
    }
    return(list(`2stpT` = s0, PcU2stpT = pct, ordr = d))
}
