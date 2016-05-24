plot.edr <-
function(x,...){

	if (!inherits(x, "edr")) 
		stop("use only with \"edr\" objects")

	ylab<-"Y"

	if (is.null(x$indices)) {
		x$indices <- x$X%*%x$matEDR
	}


	if(x$K<=2){
		for(i in 1:x$K){
			dev.new()
			xlab<-paste("Estimated index",i)
			plot(x$Y~x$indices[,i],xlab=xlab,ylab=ylab)
			lines(ksmooth(x$indices[,i],x$Y, "normal", bandwidth=0.5), col=3)	
		}
	} else {
		dev.new()
		name.var <- c("Y",paste("Estimated index",1:x$K,sep=" ")) 
		panel.smooth1<-function(x,y){panel.smooth(x,y,span=0.2)}
		pairs(cbind(x$Y,x$indices[,1:x$K]), panel=panel.smooth1, lower.panel=NULL,
				 labels=name.var)
	}

	if(x$K==2){
		X1<-x$indices[,1]
		X2<-x$indices[,2]
		X1new <- seq(min(X1), max(X1), len=50)
		X2new <- seq(min(X2), max(X2), len=50)
		newdata <- expand.grid(X1=X1new,X2=X2new)
		mod.lo<-loess(x$Y ~ X1 + X2, span=.5, degree=2)
		fit <- matrix(predict(mod.lo, newdata), 50, 50)
		persp3d(X1new, X2new, fit, xlab="Estimated index 1", ylab="Estimated index 2",
				zlab=ylab, zlim=range(c(fit,x$Y)), col="lightblue")
		points3d(X1,X2,x$Y, col = "red", pch =1)	
	}

}

