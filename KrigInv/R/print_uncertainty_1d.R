
print_uncertainty_1d <- function(model,T,type="pn",lower=0,upper=1,resolution=500,new.points=0,xlab="",ylab="",
			main="",xscale=c(0,1),show.points=TRUE,cex.main=1,cex.lab=1,
			cex.points=1,cex.axis=1,pch.points.init=17,pch.points.end=17,
			col.points.init="black",col.points.end="red",
			xaxislab=NULL,yaxislab=NULL,xaxispoint=NULL,yaxispoint=NULL,xdecal=3,ydecal=3,
      DiceViewplot=FALSE,vorobmean=FALSE){
	
	n <- model@n
	initSize <- n -new.points
	obs.X <- as.numeric(model@X)
	alpha <- NULL
	
	s <- matrix(seq(from=lower,to=upper,length=resolution),ncol=1)
	
	pred <- predict_nobias_km(object=model,newdata=s,type="UK")
	pred_obs <- predict_nobias_km(object=model,newdata=model@X,type="UK")
	
	pn <- pnorm((pred$mean - T)/pred$sd)
	pn_obs <- pnorm((pred_obs$mean - T) / pred_obs$sd)
		
	if(type=="pn"){
		myvect <- pn
		obs.Y <- pn_obs
	}else if(type=="sur"){
		myvect <- pn * (1-pn)
		obs.Y <- rep(0,times=n)
	}else if(type=="timse"){
		Wn <- 1/(sqrt(2*pi)*pred$sd)*exp(-1/2*((pred$mean-T)/pred$sd)^2)
		myvect <- (pred$sd)^2 * Wn
		obs.Y <- rep(0,times=n)
	}
	else if(type=="imse"){
		myvect <- (pred$sd)^2
		obs.Y <- rep(0,times=n)
	}else if(type=="vorob"){
    alpha <- vorob_threshold(pn=pn)
    pn_bigger_than_alpha <- (pn>alpha)+0
    pn_lower_than_alpha <- 1-pn_bigger_than_alpha
    myvect <- pn*pn_lower_than_alpha + (1-pn)*pn_bigger_than_alpha
		obs.Y <- rep(0,times=n)  
	}else{
		print("unknown value for the type argument, we take type=pn")		
		myvect <- pn
		obs.Y <- pn_obs
	}
	
	scale.x <- xscale[1] + seq(from=lower,to=upper,length=resolution) * (xscale[2]-xscale[1])
	obs.X <- xscale[1] + obs.X * (xscale[2]-xscale[1])
	axes <- is.null(xaxislab)
    
  if(DiceViewplot) par(mfrow=c(1,2))
  
	plot(x = scale.x, y = myvect, type = "l", xlab = "", ylab = "",main=main,cex.axis=cex.axis,cex.main=cex.main,cex.lab=cex.lab,axes=axes,lwd=2)
	mtext(xlab, side=1, line=xdecal,cex=cex.lab )
	mtext(ylab, side=2, line=ydecal,cex=cex.lab )
	
	if(!is.null(xaxislab)) {
		axis(side=1,at=xaxispoint,labels=xaxislab,cex.axis=cex.axis)
		axis(side=2,at=yaxispoint,labels=yaxislab,cex.axis=cex.axis)
	}
	
	if(show.points){
		points(x=obs.X[1:initSize],y=obs.Y[1:initSize], col=col.points.init, pch=pch.points.init,cex=cex.points)
	
		if (new.points!=0){
			indices <- c((initSize+1):n)
			points(x=obs.X[indices],y=obs.Y[indices], col=col.points.end, pch=pch.points.end,cex=cex.points)	
			}
	}
    
	if(vorobmean){
		if(is.null(alpha)) alpha <- vorob_threshold(pn=pn)
		lines(scale.x,rep(alpha,times=resolution),col="blue",lty=2,lwd=3)
	}
  

DiceViewplot <- FALSE #desactivated
  if(DiceViewplot){
    ymin <- min(pred$mean-3*pred$sd)
    ymax <- max(pred$mean+3*pred$sd)
    
#sectionview.km(model,ylim = c(ymin,ymax),xlim=c(lower,upper),
#title="kriging mean and variance",Xname="x",yname="f(x)",)

	  points(x=obs.X[1:initSize],y=pred_obs$mean[1:initSize], col=col.points.init, pch=pch.points.init,cex=cex.points)
	  if (new.points!=0){
    	indices <- c((initSize+1):n)
	  	points(x=obs.X[indices],y=pred_obs$mean[indices], col=col.points.end, pch=pch.points.end,cex=cex.points)	
	  }
    lines(scale.x,rep(T,times=resolution),col="black",lty=2,lwd=3)
  }
  	
	return(mean(myvect))

}
