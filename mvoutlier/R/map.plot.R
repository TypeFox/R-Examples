"map.plot" <-
function(coord, data, quan=1/2, alpha=0.025, symb=FALSE, plotmap=TRUE, map="kola.background",which.map=c(1,2,3,4),map.col=c(5,1,3,4),map.lwd=c(2,1,2,1), ... ) {

	#library(rrcov)

	# data(list=map) # is now in DESCRIPTION with LazyData: TRUE
	if(ncol(coord) != 2) stop("argument coord has to be two-dimensional")  

	rob <- covMcd(data, alpha=quan)
	dist <- mahalanobis(data, center=rob$center, cov=rob$cov)
	xarw <- arw(data, rob$center, rob$cov, alpha=alpha)
  
  	if(xarw$cn != Inf) { alpha <- sqrt(c(xarw$cn, qchisq(c(0.75,0.5,0.25),ncol(data)))) }
	else { alpha <- sqrt(qchisq(c(0.975, 0.75,0.5,0.25),ncol(data))) }
  
	if(symb==FALSE) {
		plot(coord, col=((sqrt(dist)<alpha[1])+2))
		if(plotmap == TRUE) pkb(map=map,which.map=which.map,map.col=map.col,map.lwd=map.lwd, add.plot=TRUE)
		o <- ( sqrt(dist) > min(sqrt(xarw$cn), sqrt(qchisq(0.975, dim(data)[2]) ) ) )
		l <- list(outliers = o, md=sqrt(dist))
	}
	
  
	if(symb==TRUE) {
		rd <- sqrt(dist)
		lpch <- c(3,3,16,1,1)
		lcex <- c(1.5,1,0.5,1,1.5)
		lalpha <- length(alpha)

		xs <- scale(data) - min(scale(data))
		eucl <- sqrt(apply(xs^2, 1, sum))
		rbcol <- rev(rainbow(nrow(data),start=0,end=0.7))[as.integer(cut(eucl,nrow(data),labels=1:nrow(data)))]
		
		for(j in 1:lalpha) {
			if(j==1) {
				plot(coord,type="n", ...)
				points(coord[rd>=alpha[j],],pch=lpch[j],cex=lcex[j],col=rbcol[rd>=alpha[j]])
				}
			if (j>1 & j<lalpha) points(coord[rd<alpha[j-1] & rd>=alpha[j],],cex=lcex[j],pch=lpch[j], col=rbcol[rd<alpha[j-1] & rd>=alpha[j]])
			if (j==lalpha){
           			points(coord[rd<alpha[j-1] & rd>=alpha[j],],cex=lcex[j],pch=lpch[j], col=rbcol[rd<alpha[j-1] & rd>=alpha[j]])
		        	points(coord[rd<alpha[j],],pch=lpch[j+1],cex=lcex[j+1], col=rbcol[rd<alpha[j]])
        		}
		}
		if(plotmap == TRUE) pkb(map=map,which.map=which.map,map.col=map.col,map.lwd=map.lwd, add.plot=TRUE)
		o <- ( rd > min(sqrt(xarw$cn), sqrt(qchisq(0.975, dim(data)[2]) ) ) )
		l <- list(outliers = o, md=sqrt(dist), euclidean=eucl)
	}
	l
}
