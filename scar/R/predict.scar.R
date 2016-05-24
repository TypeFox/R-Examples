predict.scar <-function(object,newdata,type = c("link", "response"),rule=1,...){
	library(stats)
	n=nrow(object$x)
	d=ncol(object$x)
        # Some checking on inputs
	if (missing(newdata)) { result = apply(object$componentfit,1,sum)+object$constant} else {
        	if(is.vector(newdata)==TRUE && d == 1) newdata=matrix(newdata,ncol=1)
		if(is.vector(newdata)==TRUE && d > 1) newdata=matrix(newdata,nrow=1)
		if(is.matrix(newdata)==FALSE||is.numeric(newdata)==FALSE) stop("Input error!")
        	if(ncol(newdata)!=d) stop("Input error: dimension mismatch.")
		total=nrow(newdata)

		fitted = matrix(0,ncol=d, nrow = total)
		for(j in 1:d){
			if (rule==1){
        			fitted[,j] = approx(object$x[,j],object$componentfit[,j],newdata[,j])$y
				xmin = min(object$x[,j])
				xmax = max(object$x[,j])
				epsilon = 1e-3
				ymin = approx(object$x[,j],object$componentfit[,j],c(xmin,xmin+epsilon),rule=2)$y
				ymax = approx(object$x[,j],object$componentfit[,j],c(xmax-epsilon,xmax),rule=2)$y
				dmin = diff(ymin)/epsilon
				dmax = diff(ymax)/epsilon
				for(i in 1:total) if(is.na(fitted[i,j]) == TRUE){
					if (newdata[i,j] < xmin) fitted[i,j] = ymin[1] + (newdata[i,j]-xmin)*dmin
					else if (newdata[i,j] > xmax) fitted[i,j] = ymax[2] + (newdata[i,j]-xmax)*dmax
				}
			} else {fitted[,j] = approx(object$x[,j],object$componentfit[,j],newdata[,j],rule=2)$y}
		}
		result = apply(fitted,1,sum)+object$constant
	}
	type <- match.arg(type)
	if (type == "link") return(result) else return (object$family$linkinv(result))
}
