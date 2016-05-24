"uni.plot" <-
function(x, symb=FALSE, quan=1/2, alpha=0.025, ...) {
	if(!is.matrix(x) && !is.data.frame(x)) stop("x must be matrix or data.frame")
	if(ncol(x) < 2) stop("x must be at least two-dimensional")
	if(ncol(x) > 10) stop("x should not be more than 10-dimensional")
	
	#library(rrcov)
	par(mfrow=c(1, ncol(x)), mai=c(0.6,0,0.6,0), oma=c(0,3,0,3))
	
	rob <- covMcd(x, alpha=quan)
	xarw <- arw(x, rob$center, rob$cov, alpha=alpha)
	
	if(xarw$cn != Inf) { alpha <- sqrt(c(xarw$cn, qchisq(c(0.75,0.5,0.25),ncol(x)))) }
	else { alpha <- sqrt(qchisq(c(0.975, 0.75,0.5,0.25),ncol(x))) }
	
	dist <- mahalanobis(x, center=rob$center, cov=rob$cov)
	sx <- matrix(NA, nrow=nrow(x), ncol=ncol(x))
	for(i in 1:ncol(x)) sx[,i] <- (x[,i]-xarw$m[i])/sqrt(xarw$c[i,i])
	r <- range(sx)
	
	if(symb == FALSE) {
		for(i in 1:ncol(x)) {
			plot(runif(nrow(x), min=-1, max=1), sx[,i], main=dimnames(x)[[2]][i], xlim=c(-1.5,1.5), ylim=c(r[1], r[2]), xlab="", ylab="Scaled Data", xaxt="n", col=(sqrt(dist)<alpha[1])+2, ...)
			par(yaxt="n")
			abline(h=0, lty="dotted")
			o <- ( sqrt(dist) > min(sqrt(xarw$cn), sqrt(qchisq(0.975, dim(x)[2]) ) ) )
			l <- list(outliers = o, md=sqrt(dist))
		}
	}
	
	if(symb == TRUE) {
		rd <- sqrt(dist)
		lpch <- c(3,3,16,1,1)
		lcex <- c(1.5,1,0.5,1,1.5)
		lalpha <- length(alpha)
		
		xs <- scale(x) - min(scale(x))
		eucl <- sqrt(apply(xs^2, 1, sum))
		rbcol <- rev(rainbow(nrow(x),start=0,end=0.7))[as.integer(cut(eucl,nrow(x),labels=1:nrow(x)))]
		
		for(i in 1:ncol(x)) {
			for(j in 1:lalpha) {
				if(j==1) {
					plot(runif(nrow(x), min=-1, max=1), sx[,i], main=dimnames(x)[[2]][i], xlim=c(-1.5,1.5), ylim=c(r[1], r[2]), xlab="", ylab="Scaled Data", xaxt="n", type="n", ...)
					par(yaxt="n")
					points(runif(nrow(x), min=-1, max=1)[rd>=alpha[j]], sx[rd>=alpha[j],i], pch=lpch[j],cex=lcex[j],col=rbcol[rd>=alpha[j]])
				}
				if (j>1 & j<lalpha) points(runif(nrow(x), min=-1, max=1)[rd<alpha[j-1] & rd>=alpha[j]], sx[rd<alpha[j-1] & rd>=alpha[j],i], cex=lcex[j],pch=lpch[j], col=rbcol[rd<alpha[j-1] & rd>=alpha[j]])
				if (j==lalpha){
           				points(runif(nrow(x), min=-1, max=1)[rd<alpha[j-1] & rd>=alpha[j]], sx[rd<alpha[j-1] & rd>=alpha[j],i], cex=lcex[j],pch=lpch[j], col=rbcol[rd<alpha[j-1] & rd>=alpha[j]])
		        		points(runif(nrow(x), min=-1, max=1)[rd<alpha[j]], sx[rd<alpha[j],i], pch=lpch[j+1],cex=lcex[j+1], col=rbcol[rd<alpha[j]])
        			}
        		}
        		abline(h=0, lty="dotted")
        	}
        o <- ( sqrt(dist) > min(sqrt(xarw$cn), sqrt(qchisq(0.975, dim(x)[2]) ) ) )
        l <- list(outliers=o, md=sqrt(dist), euclidean=eucl)
        }
        par(yaxt="s")
	l
}

