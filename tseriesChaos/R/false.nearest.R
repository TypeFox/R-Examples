#Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:15:39 +0100 (ven, 02 dic 2005) $
false.nearest <- function(series, m, d, t, rt=10, eps=sd(series)/10) {
	x <- (series - min(series)) / diff(range(series))
	eps <- eps/diff(range(series))
	res <- numeric(m)
	res2 <- numeric(m)
	for(i in 1:m) {
		a <- .C("falseNearest", series=as.double(x), length=as.integer(length(series)), m=as.integer(i), d=as.integer(d), t=as.integer(t), eps=as.double(eps), rt=as.double(rt), out=as.double(res[i]), out2=as.integer(res2[i]), PACKAGE="tseriesChaos")
		res[i] = a[["out"]]
		res2[i]= a[["out2"]]
	}
	res <- rbind(res, res2)
	res[res==(-1)] <- NA
	rownames(res) <- c("fraction", "total")
	colnames(res) <- paste("m", 1:m, sep="")
	class(res) <- "false.nearest"
	res
}

print.false.nearest <- function(x, ...) {
	print.table(x)
}

plot.false.nearest <- function(x, ...) {
	y <- x[1, !is.na(x[2,])]*100
	x <- (1:ncol(x))[!is.na(x[2,])]
	plot(x, y, type="b", main="False nearest neighbours", xlab="embedding dimension", ylab="% of false nearest neighbours", ylim=c(0,100))
}
