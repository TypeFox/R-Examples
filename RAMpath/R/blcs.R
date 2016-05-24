## bivariate latent change score model
## Johnny Zhang & Jack McArdle, 
## Created on May 26, 2012

ramLCS<-function(data, ## data to be used
	y, ## variables x to be used
	timey, ##
	ram.out=FALSE,
	betay, my0, mys, varey, vary0, varys, vary0ys,...){
	if (missing(data)) stop("No data was provided!")
	if (missing(y)) stop("Missing y variables!")
	if (!is.data.frame(data)) stop("The provided data set is not a data frame!")
	
	varname<-names(data)
	ny<-length(y)
	
	if (missing(timey)) timey<-1:ny
	
	My<-max(timey)
	
	## for Y	
	model<-'y0 =~ 1*y1 \n '
	
	for (i in 2:My){
		model<-paste(model, "y",i,"~1*y",(i-1),"\n ", sep="")
	}
	
	for (i in 2:My){
		model<-paste(model, "dy",i,"=~1*y",i,"\n ", sep="")
	}
	
	for (i in 2:My){
		model<-paste(model, "dy",i,"~betay*y", (i-1), sep="")
		if (!missing(betay)){
			model<-paste(model, "+start(", betay, ")*y", (i-1), "\n ", sep="")
		}else{
			model<-paste(model, "\n ")
		}
	}
	for (i in 2:My){
		model<-paste(model, "ys=~1*dy", i, "\n ", sep="")
	}
	
	for (i in 2:My){
		model<-paste(model, "dy",i,"~~0*dy",i,"\n ", sep="")
	}
	
	for (i in 1:My){
		model<-paste(model, "y",i,"~~0*y",i,"\n ", sep="")
	}
	
	model<-paste(model, "ys~~vary0ys*y0")
	if (!missing(vary0ys)){
		model<-paste(model, "+start(",vary0ys,")*y0\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	model<-paste(model, "y0~~vary0*y0")
	if (!missing(vary0)){
		model<-paste(model, "+start(",vary0,")*y0\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	model<-paste(model, "ys~~varys*ys")
	if (!missing(varys)){
		model<-paste(model, "+start(",varys,")*ys\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	model<-paste(model, "ys~mys*1")
	if (!missing(mys)){
		model<-paste(model, "+start(",mys,")*1\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	model<-paste(model, "y0~my0*1")
	if (!missing(my0)){
		model<-paste(model, "+start(",my0,")*1\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	for (i in 1:My){
		model<-paste(model, "y",i,"~0*1\n ", sep="")
	}
	for (i in 2:My){
		model<-paste(model, "dy",i,"~0*1\n ", sep="")
	}
	
	
	## for observed data part
	if (ny < My){
		ally<-1:My
		missy<-ally[-timey]
		for (i in 1:ny){
			model<-paste(model, "y",timey[i],"=~1*", varname[y[i]],"\n", sep="")				
		}
		for (j in missy){
			model<-paste(model, "y",j,"=~0*",varname[y[1]],"\n", sep="")	
		}
		for (i in 1:ny){
			model<-paste(model, varname[y[i]], "~0*1\n", sep="")		
		}
		for (i in 1:ny){
			model<-paste(model, varname[y[i]], "~~varey*",varname[y[i]], sep="")
			if (!missing(varey)){
				model<-paste(model, "+start(",varey,")*",varname[y[i]], "\n ", sep="")
			}else{
				model<-paste(model, "\n ")
			}		
		}
	}else{
		for (i in 1:ny){
			model<-paste(model, "y",timey[i],"=~1*", varname[y[i]],"\n", sep="")				
		}

		for (i in 1:ny){
			model<-paste(model, varname[y[i]], "~0*1\n", sep="")		
		}
		for (i in 1:ny){
			model<-paste(model, varname[y[i]], "~~varey*",varname[y[i]], sep="")
			if (!missing(varey)){
				model<-paste(model, "+start(",varey,")*",varname[y[i]], "\n ", sep="")
			}else{
				model<-paste(model, "\n ")
			}		
		}
	}
	fitModel<-lavaan(model=model, data=data, ...)
	summary(fitModel, fit.measures=TRUE)
	if (ram.out){ 
		ram=lavaan2ram(fitModel)
		invisible((list(model=model, lavaan=fitModel, ram=ram, data=data)))
	}else{
	  ram=lavaan2ram(fitModel,ram.out=FALSE)
		invisible((list(model=model, lavaan=fitModel, ram=ram, data=data)))
	}
}


ramBLCS<-function(
	data, ## data to be used
	y, ## variables x to be used
	x, ## variables y to be used
	timey, ##
	timex, ##
	ram.out=FALSE,
	betax, betay, gammax, gammay, mx0, mxs, my0, mys, varex, varey, varx0, vary0, varxs, varys, varx0y0, varx0xs, vary0ys, varx0ys, vary0xs, varxsys,...){
	if (missing(data)) stop("No data was provided!")
	if (missing(y)) stop("Missing y variables!")
	if (!is.data.frame(data)) stop("The provided data set is not a data frame!")
	
	varname<-names(data)
	ny<-length(y)
	
	if (missing(timey)) timey<-1:ny
	
	My<-max(timey)
	
	## for Y	
	model<-'y0 =~ 1*y1 \n '
	
	for (i in 2:My){
		model<-paste(model, "y",i,"~1*y",(i-1),"\n ", sep="")
	}
	
	for (i in 2:My){
		model<-paste(model, "dy",i,"=~1*y",i,"\n ", sep="")
	}
	
	for (i in 2:My){
		model<-paste(model, "dy",i,"~betay*y", (i-1), sep="")
		if (!missing(betay)){
			model<-paste(model, "+start(", betay, ")*y", (i-1), "\n ", sep="")
		}else{
			model<-paste(model, "\n ")
		}
	}
	for (i in 2:My){
		model<-paste(model, "ys=~1*dy", i, "\n ", sep="")
	}
	
	for (i in 2:My){
		model<-paste(model, "dy",i,"~~0*dy",i,"\n ", sep="")
	}
	
	for (i in 1:My){
		model<-paste(model, "y",i,"~~0*y",i,"\n ", sep="")
	}
	
	model<-paste(model, "ys~~vary0ys*y0")
	if (!missing(vary0ys)){
		model<-paste(model, "+start(",vary0ys,")*y0\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	model<-paste(model, "y0~~vary0*y0")
	if (!missing(vary0)){
		model<-paste(model, "+start(",vary0,")*y0\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	model<-paste(model, "ys~~varys*ys")
	if (!missing(varys)){
		model<-paste(model, "+start(",varys,")*ys\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	model<-paste(model, "ys~mys*1")
	if (!missing(mys)){
		model<-paste(model, "+start(",mys,")*1\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	model<-paste(model, "y0~my0*1")
	if (!missing(my0)){
		model<-paste(model, "+start(",my0,")*1\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	for (i in 1:My){
		model<-paste(model, "y",i,"~0*1\n ", sep="")
	}
	for (i in 2:My){
		model<-paste(model, "dy",i,"~0*1\n ", sep="")
	}
	
	## for observed data part
	if (ny < My){
		ally<-1:My
		missy<-ally[-timey]
		for (i in 1:ny){
			model<-paste(model, "y",timey[i],"=~1*", varname[y[i]],"\n", sep="")				
		}
		for (j in missy){
			model<-paste(model, "y",j,"=~0*",varname[y[1]],"\n", sep="")	
		}
		for (i in 1:ny){
			model<-paste(model, varname[y[i]], "~0*1\n", sep="")		
		}
		for (i in 1:ny){
			model<-paste(model, varname[y[i]], "~~varey*",varname[y[i]], sep="")
			if (!missing(varey)){
				model<-paste(model, "+start(",varey,")*",varname[y[i]], "\n ", sep="")
			}else{
				model<-paste(model, "\n ")
			}		
		}
	}else{
		for (i in 1:ny){
			model<-paste(model, "y",timey[i],"=~1*", varname[y[i]],"\n", sep="")				
		}

		for (i in 1:ny){
			model<-paste(model, varname[y[i]], "~0*1\n", sep="")		
		}
		for (i in 1:ny){
			model<-paste(model, varname[y[i]], "~~varey*",varname[y[i]], sep="")
			if (!missing(varey)){
				model<-paste(model, "+start(",varey,")*",varname[y[i]], "\n ", sep="")
			}else{
				model<-paste(model, "\n ")
			}		
		}
	}
	
	if (missing(x)) stop("Missing x variables!")

	nx<-length(x)
	
	if (missing(timex)) timex<-1:nx
	
	Mx<-max(timex)
	
	## for x	
	model<-paste(model,'x0 =~ 1*x1 \n ')
	
	for (i in 2:Mx){
		model<-paste(model, "x",i,"~1*x",(i-1),"\n ", sep="")
	}
	
	for (i in 2:Mx){
		model<-paste(model, "dx",i,"=~1*x",i,"\n ", sep="")
	}
	
	for (i in 2:Mx){
		model<-paste(model, "dx",i,"~betax*x", (i-1), sep="")
		if (!missing(betax)){
			model<-paste(model, "+start(", betax, ")*x", (i-1), "\n ", sep="")
		}else{
			model<-paste(model, "\n ")
		}
	}
	for (i in 2:Mx){
		model<-paste(model, "xs=~1*dx", i, "\n ", sep="")
	}
	
	for (i in 2:Mx){
		model<-paste(model, "dx",i,"~~0*dx",i,"\n ", sep="")
	}
	
	for (i in 1:Mx){
		model<-paste(model, "x",i,"~~0*x",i,"\n ", sep="")
	}
	
	model<-paste(model, "xs~~varx0xs*x0")
	if (!missing(varx0xs)){
		model<-paste(model, "+start(",varx0xs,")*x0\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	model<-paste(model, "x0~~varx0*x0")
	if (!missing(varx0)){
		model<-paste(model, "+start(",varx0,")*x0\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	model<-paste(model, "xs~~varxs*xs")
	if (!missing(varxs)){
		model<-paste(model, "+start(",varxs,")*xs\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	model<-paste(model, "xs~mxs*1")
	if (!missing(mxs)){
		model<-paste(model, "+start(",mxs,")*1\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	model<-paste(model, "x0~mx0*1")
	if (!missing(mx0)){
		model<-paste(model, "+start(",mx0,")*1\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	for (i in 1:Mx){
		model<-paste(model, "x",i,"~0*1\n ", sep="")
	}
	for (i in 2:Mx){
		model<-paste(model, "dx",i,"~0*1\n ", sep="")
	}
	
	## for observed data part
	if (nx < Mx){
		allx<-1:Mx
		missx<-allx[-timex]
		for (i in 1:nx){
			model<-paste(model, "x",timex[i],"=~1*", varname[x[i]],"\n", sep="")				
		}
		for (j in missx){
			model<-paste(model, "x",j,"=~0*",varname[x[1]],"\n", sep="")	
		}
		for (i in 1:nx){
			model<-paste(model, varname[x[i]], "~0*1\n", sep="")		
		}
		for (i in 1:nx){
			model<-paste(model, varname[x[i]], "~~varex*",varname[x[i]], sep="")
			if (!missing(varex)){
				model<-paste(model, "+start(",varex,")*",varname[x[i]], "\n ", sep="")
			}else{
				model<-paste(model, "\n ")
			}		
		}
	}else{
		for (i in 1:nx){
			model<-paste(model, "x",timex[i],"=~1*", varname[x[i]],"\n", sep="")				
		}

		for (i in 1:nx){
			model<-paste(model, varname[x[i]], "~0*1\n", sep="")		
		}
		for (i in 1:nx){
			model<-paste(model, varname[x[i]], "~~varex*",varname[x[i]], sep="")
			if (!missing(varex)){
				model<-paste(model, "+start(",varex,")*",varname[x[i]], "\n ", sep="")
			}else{
				model<-paste(model, "\n ")
			}		
		}
	}
	
	## coupling effects
	for (i in 2:My){
		model<-paste(model, "dy",i,"~gammax*x",i-1, sep="")
		if (!missing(gammax)){
			model<-paste(model, "+start(",gammax,")*x",i-1, "\n ", sep="")
		}else{
			model<-paste(model, "\n ")
		}
	}
	
	for (i in 2:Mx){
		model<-paste(model, "dx",i,"~gammay*y",i-1, sep="")
		if (!missing(gammay)){
			model<-paste(model, "+start(",gammay,")*y",i-1, "\n ", sep="")
		}else{
			model<-paste(model, "\n ")
		}
	}
	
	model<-paste(model, "x0~~varx0y0*y0")
	if (!missing(varx0y0)){
		model<-paste(model, "+start(",varx0y0,")*y0\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	model<-paste(model, "x0~~varx0ys*ys")
	if (!missing(varx0ys)){
		model<-paste(model, "+start(",varx0ys,")*ys\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	model<-paste(model, "y0~~vary0xs*xs")
	if (!missing(vary0xs)){
		model<-paste(model, "+start(",vary0xs,")*xs\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	model<-paste(model, "xs~~varxsys*ys")
	if (!missing(varxsys)){
		model<-paste(model, "+start(",varxsys,")*ys\n ")
	}else{
		model<-paste(model, "\n ")
	}
	
	fitModel<-lavaan(model=model, data=data, ...)

	summary(fitModel, fit.measures=TRUE)
	if (ram.out){ 
		ram=lavaan2ram(fitModel)
		blcs<-list(model=model, lavaan=fitModel, ram=ram, info=list(y=y,x=x), data=data)
		class(blcs)<-'blcs'
		invisible(blcs)
	}else{
		ram=lavaan2ram(fitModel,ram.out=FALSE)
		blcs<-list(model=model, lavaan=fitModel, ram=ram, info=list(y=y,x=x), data=data)
		class(blcs)<-'blcs'
		invisible(blcs)
	}
}


ramVF<-function(ramout, ylim, xlim, ninterval=10, scale=.1, length=.25, scatter=TRUE, n=20, alpha=.95, ...){
	ind<-which(ramout$lavaan@ParTable$label=="betax")[1]
	betax<-ramout$lavaan@Fit@est[ind]
	
	ind<-which(ramout$lavaan@ParTable$label=="gammax")[1]
	gammax<-ramout$lavaan@Fit@est[ind]
	
	ind<-which(ramout$lavaan@ParTable$label=="betay")[1]
	betay<-ramout$lavaan@Fit@est[ind]
	
	ind<-which(ramout$lavaan@ParTable$label=="gammay")[1]
	gammay<-ramout$lavaan@Fit@est[ind]
	
	ind<-which(ramout$lavaan@ParTable$label=="mxs")[1]
	mux<-ramout$lavaan@Fit@est[ind]
	
	ind<-which(ramout$lavaan@ParTable$label=="mys")[1]
	muy<-ramout$lavaan@Fit@est[ind]
	
	x<-seq(xlim[1],xlim[2], length=ninterval)
	y<-seq(ylim[1],ylim[2], length=ninterval)
	xy<-expand.grid(x,y)
	
	x<-xy[,1]
	y<-xy[,2]
	
	dx<-mux + betax*x + gammay*y
	dy<-muy + betay*y + gammax*x
	
	x1<-x+scale*dx
	y1<-y+scale*dy
	
	plot(x,y,type='n', ...)
	arrows(x,y,x1,y1,length=length,...)	
	
	if (scatter){
		alldata<-ramout$lavaan@Data@X[[1]]
		xdata<-c(alldata[1:n, ramout$info$x])
		ydata<-c(alldata[1:n, ramout$info$y])
		points(xdata, ydata, col='grey', ...)
		
		## add confidence interval
		##
		xall <- c(alldata[,ramout$info$x])
		yall <- c(alldata[,ramout$info$y])
			       
		       lines(ellipse(cor(xall, yall), level=alpha, scale=c(sd(xall),sd(yall)), centre=c(mean(xall),mean(yall))), lwd=1.5, col="green")
	}
	
	invisible(cbind(x,y,x1,y1))
}

plot.blcs<-function(x, ylim, xlim, ninterval=10, scale=.1, length=.25, scatter=TRUE, n=20, alpha=.95, ...){
	ramout<-x
	ind<-which(ramout$lavaan@ParTable$label=="betax")[1]
	betax<-ramout$lavaan@Fit@est[ind]
	
	ind<-which(ramout$lavaan@ParTable$label=="gammax")[1]
	gammax<-ramout$lavaan@Fit@est[ind]
	
	ind<-which(ramout$lavaan@ParTable$label=="betay")[1]
	betay<-ramout$lavaan@Fit@est[ind]
	
	ind<-which(ramout$lavaan@ParTable$label=="gammay")[1]
	gammay<-ramout$lavaan@Fit@est[ind]
	
	ind<-which(ramout$lavaan@ParTable$label=="mxs")[1]
	mux<-ramout$lavaan@Fit@est[ind]
	
	ind<-which(ramout$lavaan@ParTable$label=="mys")[1]
	muy<-ramout$lavaan@Fit@est[ind]
	
	x<-seq(xlim[1],xlim[2], length=ninterval)
	y<-seq(ylim[1],ylim[2], length=ninterval)
	xy<-expand.grid(x,y)
	
	x<-xy[,1]
	y<-xy[,2]
	
	dx<-mux + betax*x + gammay*y
	dy<-muy + betay*y + gammax*x
	
	x1<-x+scale*dx
	y1<-y+scale*dy
	
	plot(x,y,type='n', ...)
	arrows(x,y,x1,y1,length=length,...)	
	
	if (scatter){
		alldata<-ramout$lavaan@Data@X[[1]]
		xdata<-c(alldata[1:n, ramout$info$x])
		ydata<-c(alldata[1:n, ramout$info$y])
		points(xdata, ydata, col='grey', ...)
		
		## add confidence interval
		##
		xall <- c(alldata[,ramout$info$x])
		yall <- c(alldata[,ramout$info$y])
				       
		lines(ellipse(cor(xall, yall), level=alpha, scale=c(sd(xall),sd(yall)), centre=c(mean(xall),mean(yall))), lwd=1.5, col="green")
	}
	
	invisible(cbind(x,y,x1,y1))
}
