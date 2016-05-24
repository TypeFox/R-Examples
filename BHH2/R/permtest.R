"permtest" <-
function (x, y = NULL) 
{
    if (is.null(y)) {
        mx <- mean(x)
        n <- length(x)
        t.obs <- mx/sqrt(var(x)/n)
        N <- 2^n
        mat <- matrix(NA, nrow = N, ncol = n)
        for (j in 1:n) {
            m <- 2^j
            mat[, j] <- rep(c(rep(-1, N/m), rep(+1, N/m)), m/2)
        }
        d <- as.numeric(mat %*% x/n)
        k <- length(which(d > mx))
        return(c(N = N, t.obs = t.obs, "t-Dist-P(>t)" = 1 - pt(t.obs, 
            n - 1), "PermDist-P(>t)" = k/N))
    }
    else {
        x <- x[!is.na(x)]
        y <- y[!is.na(y)]
        nx <- length(x)
        ny <- length(y)
        S2x <- sum(x^2) - sum(x)^2/nx
        S2y <- sum(y^2) - sum(y)^2/ny
        t.stat <- function(x, y) {
            (mean(x) - mean(y))/sqrt((S2x + S2y)/(nx + ny - 2) * 
                (1/nx + 1/ny))
        }
        f.stat <- function(x, y) {
            (S2x/(nx - 1))/(S2y/(ny - 1))
        }
        t.obs <- t.stat(x, y)
        f.obs <- f.stat(x, y)
        z <- c(x, y)
        n <- nx + ny
        mat <- subsets(n, nx)
        N <- nrow(mat)
        kt <- kf <- 0
        for (i in 1:nrow(mat)) {
            x <- z[mat[i, ]]
            y <- z[-mat[i, ]]
            S2x <- sum(x^2) - sum(x)^2/nx
            S2y <- sum(y^2) - sum(y)^2/ny
            if (t.obs < t.stat(x, y)) 
                kt <- kt + 1
            if (f.obs < f.stat(x, y)) 
                kf <- kf + 1
        }
        return(c(N = N, t.obs = t.obs, "t-Dist:P(>t)" = 1 - pt(t.obs, 
            nx + ny - 2), "PermDist:P(>t)" = kt/N, F.obs = f.obs, 
            "F-Dist:P(>F)" = 1 - pf(f.obs, nx - 1, ny - 1), "PermDist:P(>F)" = kf/N))
    }
}
