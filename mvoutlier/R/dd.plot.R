"dd.plot" <-
function(x, quan=1/2, alpha=0.025, ...) {
	if(!is.matrix(x) && !is.data.frame(x)) stop("x must be matrix or data.frame")
	#library(rrcov)
	rob <- covMcd(x, alpha=quan)
	xarw <- arw(x, rob$center, rob$cov, alpha=alpha)
	
	distcla <- sqrt(mahalanobis(x, center=apply(x, 2, mean), cov=cov(x)))
	distrob <- sqrt(mahalanobis(x, center=rob$center, cov=rob$cov))
	plot(distcla, distrob, main="Distance-Distance Plot", xlab="Mahalanobis Distance", ylab="Robust Distance", type="n", ...)
	
	
	if(xarw$cn != Inf) { alpha <- sqrt(c(xarw$cn, qchisq(c(0.75,0.5,0.25),ncol(x)))) }
	else { alpha <- sqrt(qchisq(c(0.975, 0.75,0.5,0.25),ncol(x))) }
	abline(h=alpha[1])
	abline(v=alpha[1])
	abline(a=0, b=1)
	lpch <- c(3,3,16,1,1)
	lcex <- c(1.5,1,0.5,1,1.5)
	lalpha <- length(alpha)
	
	xs <- scale(x) - min(scale(x))
	eucl <- sqrt(apply(xs^2, 1, sum))
	rbcol <- rev(rainbow(nrow(x),start=0,end=0.7))[as.integer(cut(eucl,nrow(x),labels=1:nrow(x)))]
	rd <- distrob
	
	for(j in 1:lalpha) {
		if(j==1) {
			points(distcla[rd>=alpha[j]], distrob[rd>=alpha[j]], pch=lpch[j],cex=lcex[j],col=rbcol[rd>=alpha[j]])
		}
		if (j>1 & j<lalpha) points(distcla[rd<alpha[j-1] & rd>=alpha[j]], distrob[rd<alpha[j-1] & rd>=alpha[j]], cex=lcex[j],pch=lpch[j], col=rbcol[rd<alpha[j-1] & rd>=alpha[j]])
		if (j==lalpha){
        		points(distcla[rd<alpha[j-1] & rd>=alpha[j]], distrob[rd<alpha[j-1] & rd>=alpha[j]], cex=lcex[j],pch=lpch[j], col=rbcol[rd<alpha[j-1] & rd>=alpha[j]])
		       	points(distcla[rd<alpha[j]], distrob[rd<alpha[j]], pch=lpch[j+1],cex=lcex[j+1], col=rbcol[rd<alpha[j]])
        	}
        }
        o <- ( rd > min(sqrt(xarw$cn), sqrt(qchisq(0.975, dim(x)[2]) ) ) )
        l <- list(outliers = o, md.cla = distcla, md.rob=distrob)
        l
}

