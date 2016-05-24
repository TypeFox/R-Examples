nptol.int <- function (x, alpha = 0.05, P = 0.99, side = 1, method = c("WILKS", 
    "WALD", "HM", "YM"), upper = NULL, lower = NULL) 
{
    n <- length(x)
    x.sort <- sort(x)
    method <- match.arg(method)
    if (is.null(upper)) 
        upper <- max(x)
    if (is.null(lower)) 
        lower <- min(x)
    if (method == "WILKS") {
        if (side == 2) {
            if (floor((n + 1)/2) == ((n + 1)/2)) 
                up <- ((n + 1)/2) - 1
            else up <- floor((n + 1)/2)
            r <- 1:up
            out2 <- pbeta(P, n - 2 * r + 1, 2 * r, lower.tail = FALSE) - 
                (1 - alpha)
            ti2 <- cbind(r, out2)
            temp2 <- matrix(ti2[(ti2[, 2] > 0), ], ncol = 2)
            if (nrow(temp2) == 0) {
                lower <- lower
                upper <- upper
            }
            else {
                mins2 <- min(temp2[, 2])
                temp2 <- matrix(temp2[temp2[, 2] == mins2, ], 
                  ncol = 2)
                r = temp2[, 1]
                lower <- x.sort[r]
                upper <- x.sort[n - r + 1]
            }
        }
        if (side == 1) {
            r <- qbinom(alpha, size = n, prob = 1 - P)
            s <- n - r + 1
            if (r < 1) {
                lower <- lower
            }
            else lower <- x.sort[r]
            if (s > n) {
                upper <- upper
            }
            else upper <- x.sort[s]
            temp <- data.frame(cbind(alpha, P, lower, upper))
            colnames(temp) <- c("alpha", "P", "1-sided.lower", 
                "1-sided.upper")
        }
        else {
            temp <- data.frame(cbind(alpha, P, lower, upper))
            colnames(temp) = c("alpha", "P", "2-sided.lower", 
                "2-sided.upper")
        }
    }
    if (method == "WALD") {
        t <- NULL
        s <- NULL
        for (i in 2:n) {
            s <- c(s, 1:(i - 1))
            t <- c(t, rep(i, i - 1))
        }
        if (side == 1) {
            r <- qbinom(alpha, size = n, prob = 1 - P)
            s <- n - r + 1
            if (r < 1) {
                lower <- lower
            }
            else lower <- x.sort[r]
            if (s > n) {
                upper <- upper
            }
            else upper <- x.sort[s]
            temp <- data.frame(cbind(alpha, P, lower, upper))
        }
        else {
            out3 <- pbeta(P, t - s, n - t + s + 1, lower.tail = FALSE) - 
                (1 - alpha)
            ti3 <- cbind(s, t, out3)
            temp3 <- matrix(ti3[(ti3[, 3] > 0), ], ncol = 3)
            if (nrow(temp3) == 0) {
                lower <- lower
                upper <- upper
            }
            else {
                mins3 <- min(temp3[, 3])
                out5 <- matrix(temp3[temp3[, 3] == mins3, ], ncol = 3)
                s <- out5[, 1]
                t <- out5[, 2]
                lower <- x.sort[s]
                upper <- x.sort[t]
            }
        }
        if (side == 1) {
            temp <- data.frame(cbind(alpha, P, lower, upper))
            colnames(temp) <- c("alpha", "P", "1-sided.lower", 
                "1-sided.upper")
        }
        else {
            temp <- data.frame(cbind(alpha, P, lower, upper))
            colnames(temp) <- c("alpha", "P", "2-sided.lower", 
                "2-sided.upper")
        }
    }
    if (method == "HM") {
        ind <- 0:n
        out <- pbinom(ind, n, P) - (1 - alpha)
        ti <- cbind(ind, out)
        temp <- matrix(ti[(ti[, 2] > 0), ], ncol = 2)
        mins <- min(temp[, 2])
        HM.ind <- temp[temp[, 2] == mins, 1]
        diff <- n - HM.ind
        if (side == 2) {
            if (diff == 0 | floor(diff/2) == 0) {
                if (lower) 
                  x.sort <- c(lower, x.sort)
                if (upper) 
                  x.sort <- c(x.sort, upper)
                HM = cbind(1, length(x.sort))
            }
            else {
                if (floor(diff/2) == diff/2) {
                  v1 <- v2 <- diff/2
                }
                else {
                  v1 <- c(floor(diff/2), ceiling(diff/2))
                  v2 <- v1 + c(1, -1)
                }
                HM <- cbind(v1, n - v2 + 1)
            }
            temp <- data.frame(cbind(alpha, P, x.sort[HM[, 1]], 
                x.sort[HM[, 2]]))
	if (sum(dim(HM) == c(2, 2))==2){
        	if ((x.sort[HM[1, 1]] == x.sort[HM[2, 1]]) & (x.sort[HM[1, 2]] == x.sort[HM[2, 2]])) temp <- temp[1, ]
	}
            colnames(temp) <- c("alpha", "P", "2-sided.lower", 
                "2-sided.upper")
        }
        else {
            l <- 0:n
            u <- 1:(n + 1)
            low.temp <- cbind(l, ((1 - pbinom(l - 1, n, 1 - P)) - 
                (1 - alpha)))
            l <- matrix(low.temp[low.temp[, 2] > 0, ], ncol = 2)
            l <- l[which.max(l[, 1]), ][1]
            if (l > 0) 
                lower <- x.sort[l]
            up.temp <- cbind(u, (pbinom(u - 1, n, P)) - (1 - alpha))
            u <- matrix(up.temp[up.temp[, 2] > 0, ], ncol = 2)
            u <- u[which.min(u[, 1]), ][1]
            if (u < (n + 1)) 
                upper <- x.sort[u]
            temp <- data.frame(cbind(alpha, P, lower, upper))
            colnames(temp) <- c("alpha", "P", "1-sided.lower", 
                "1-sided.upper")
        }
    }
    if (method == "YM") {
		n.min <- as.numeric(distfree.est(alpha = alpha, P = P, side = side))
		if(side == 1){
			if(n<n.min) temp <- extrap(x=x, alpha=alpha, P=P) else temp <- interp(x=x, alpha=alpha, P=P)
		} else{
			temp <- two.sided(x=x, alpha=alpha, P=P)		
		}
    }
	
    temp
}


