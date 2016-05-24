num.sv2 <-
function (dat, mod, method = c("be", "leek"), vfilter = NULL, 
    B = 20, sv.sig = 0.1, seed = NULL) 
{
modefunc <- function(x){
	return(as.numeric(names(sort(-table(x)))[1]))
}

    if (!is.null(vfilter)) {
        if (vfilter < 100 | vfilter > dim(dat)[1]) {
            stop(paste("The number of genes used in the analysis must be between 100 and", 
                dim(dat)[1], "\n"))
        }
        tmpv = rowVars(dat)
        ind = which(rank(-tmpv) < vfilter)
        dat = dat[ind, ]
    }
    method <- match.arg(method)
    if (method == "be") {
        if (!is.null(seed)) {
            set.seed(seed)
        }
        warn <- NULL
        n <- ncol(dat)
        m <- nrow(dat)
        H <- mod %*% solve(t(mod) %*% mod) %*% t(mod)
        res <- dat - t(H %*% t(dat))
        uu <- fast.svd(res)
        ndf <- n - ceiling(sum(diag(H)))
        dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
        dstat0 <- matrix(0, nrow = B, ncol = ndf)
        for (i in 1:B) {
            res0 <- t(apply(res, 1, sample, replace = FALSE))
            res0 <- res0 - t(H %*% t(res0))
            uu0 <- fast.svd(res0)
            dstat0[i, ] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
        }
        psv <- rep(1, n)
        for (i in 1:ndf) {
            psv[i] <- mean(dstat0[, i] >= dstat[i])
        }
        for (i in 2:ndf) {
            psv[i] <- max(psv[(i - 1)], psv[i])
        }
        nsv <- sum(psv <= sv.sig)
        return(as.numeric(list(n.sv = nsv)))
    }
    else {
        dat <- as.matrix(dat)
        dims <- dim(dat)
        a <- seq(0, 2, length = 100)
        n <- floor(dims[1]/10)
        rhat <- matrix(0, nrow = 100, ncol = 10)
        P <- (diag(dims[2]) - mod %*% solve(t(mod) %*% mod) %*% 
            t(mod))
        for (j in 1:10) {
            dats <- dat[1:(j * n), ]
            ee <- eigen(t(dats) %*% dats)
            sigbar <- ee$values[dims[2]]/(j * n)
            R <- dats %*% P
            wm <- (1/(j * n)) * t(R) %*% R - P * sigbar
            ee <- eigen(wm)
            v <- c(rep(T, 100), rep(F, dims[2]))
            v <- v[order(c(a * (j * n)^(-1/3) * dims[2], ee$values), 
                decreasing = TRUE)]
            u <- 1:length(v)
            w <- 1:100
            rhat[, j] <- rev((u[v == TRUE] - w))
        }
        ss <- rowVars(rhat)
        bumpstart <- which.max(ss > (2 * ss[1]))
        start <- which.max(c(rep(1e+05, bumpstart), ss[(bumpstart + 
            1):100]) < 0.5 * ss[1])
        finish <- which.max(ss * c(rep(0, start), rep(1, 100 - 
            start)) > ss[1])
        if (finish == 1) {
            finish <- 100
        }
        n.sv <- modefunc(rhat[start:finish, 10])
        return(n.sv)
        print(method)
    }
}
