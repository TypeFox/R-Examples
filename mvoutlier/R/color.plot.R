"color.plot" <-
function(x, quan=1/2, alpha=0.025, ...)  {

	#library(rrcov)
	if(!is.matrix(x) && !is.data.frame(x)) stop("x has to be matrix or data.frame")
	if(ncol(x) != 2) stop("x has to be two-dimensional")

	n <- nrow(x)
	rob <- covMcd(x, alpha=quan)
	xarw <- arw(x, rob$center, rob$cov, alpha=alpha)
	xs <- scale(x) - min(scale(x))
	eucl <- sqrt(xs[,1]^2 + xs[,2]^2)
	rbcol <- rev(rainbow(nrow(x),start=0,end=0.7))[as.integer(cut(eucl,nrow(x),labels=1:nrow(x)))]

	covr <- rob$cov
	mer <- rob$center

	covr.svd <- svd(covr, nv = 0)
	rr <- covr.svd[["u"]] %*% diag(sqrt(covr.svd[["d"]]))

	m <- 1000
	if(xarw$cn != Inf) { alpha <- sqrt(c(xarw$cn, qchisq(c(0.75,0.5,0.25),ncol(x)))) }
	else { alpha <- sqrt(qchisq(c(0.975, 0.75,0.5,0.25),ncol(x))) }
	lpch <- c(3,3,16,1,1)
	lcex <- c(1.5,1,0.5,1,1.5)
	lalpha <- length(alpha)

	rd <- sqrt(mahalanobis(x, mer, covr))

	for(j in 1:lalpha) {
        	e1 <- cos(c(0:m)/m * 2 * pi) * alpha[j]
	        e2 <- sin(c(0:m)/m * 2 * pi) * alpha[j]
        	e <- cbind(e1, e2)
        	ttr <- t(rr %*% t(e)) + rep(1, m + 1) %o% mer
	        if(j == 1) {
        	        xmax <- max(c(x[, 1],ttr[,1]))
        	        xmin <- min(c(x[, 1],ttr[,1]))
			ymax <- max(c(x[, 2],ttr[,2]))
	                ymin <- min(c(x[, 2],ttr[,2]))
        	        plot(x, xlab = "x", ylab = "y", xlim = c(xmin, xmax), ylim = c(ymin, ymax),type="n", main="Color according to Euclidean distance", ...)
			points(x[rd>=alpha[j],],pch=lpch[j],cex=lcex[j],col=rbcol[rd>=alpha[j]])
	        }
       	        if (j>1 & j<lalpha) points(x[rd<alpha[j-1] & rd>=alpha[j],],cex=lcex[j],pch=lpch[j], col=rbcol[rd<alpha[j-1] & rd>=alpha[j]])
        	if (j==lalpha){
           		points(x[rd<alpha[j-1] & rd>=alpha[j],],cex=lcex[j],pch=lpch[j], col=rbcol[rd<alpha[j-1] & rd>=alpha[j]])
		        points(x[rd<alpha[j],],pch=lpch[j+1],cex=lcex[j+1], col=rbcol[rd<alpha[j]])
        	}
 		lines(ttr[, 1], ttr[, 2], lty=3)
	}
	
	o <- ( rd > min(sqrt(xarw$cn), sqrt(qchisq(0.975, dim(x)[2]) ) ) )
	l <- list(outliers = o, md = rd, euclidean = eucl)
	l
}

