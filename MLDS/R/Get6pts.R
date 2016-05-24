`Get6pts` <-
function(x, nrep, ...) {
# x object of class mlds
	wh <- function(xx) {
		w <- which(xx == b1)
		if (length(w) > 0) w else 
			as.integer(rep(NA, nrep))
	}
	n <- length(x$pscale)
	ii <- expand.grid(seq(dim(fs <-  t(combn(seq(n), 3)))[1]),
				seq(dim(ss <-  t(combn(seq(2, n), 3)))[1]))
	fl <- cbind(fs[ii$Var1,], ss[ii$Var2,])
	fc <- (fl[, 1] < fl[, 4])
	fl <- fl[fc, ] # all 6-tuples w/ a1 < a2

	if (x$method == "glm") {
		dd <- as.mlds.df(x$obj$data)
		} else
	{ dd <- x$data }
	
	a <- fl[, c(1, 2, 4, 5)]
	b <- fl[, c(2, 3, 5, 6)]
	e <- fl[, c(1, 3, 4, 6)] 
	
	base <- 10^seq(3, 0)
	b1 <- as.matrix((dd[, -1]) - 1) %*% base #codage of 4-tuples
	A <- as.vector(sapply(drop((a - 1) %*% base), wh))
	B <- as.vector(sapply(drop((b - 1) %*% base), wh))
	E <- as.vector(sapply(drop((e - 1) %*% base), wh))
	cc <- data.frame(A = A, B = B, E = E)
	cc <- cc[complete.cases(cc), ] #keep only if all 3 present
	Six.Pts <- lapply(cc, function(x, y) y[unlist(x), ] , y = dd)
	attr(Six.Pts, "indices") <- cc
    Six.Pts
}

