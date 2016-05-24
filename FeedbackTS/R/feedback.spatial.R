 # # # # # # # # # # # # # # # # # #
# # #     spatial feedback          # # #
 # # #         functions            # # #
  # # # # # # # # # # # # # # # # # 


map.statistic=function(coordinates,statistic,region,cex.circles=c(3,0.2),legend,
	main=NULL,add=FALSE){
  	factmult=1/max(abs(statistic))*cex.circles[1]
  	factaddi=cex.circles[2]
  	z=statistic*factmult
	map('worldHires',region$border,xlim=region$xlim,ylim=region$ylim,col=1,add=add)
	if(length(main)>0){
		mtext(main)
	}
  	points(coordinates,cex=factaddi+abs(z),col=1+(z<0)+3*(z>0),lwd=1)
  	if(length(legend$x)>0){
  		points(legend$x,legend$y,col=c(1,1,1,4,2),
        cex=factaddi+round(c(0,max(abs(z))/2,max(abs(z)),rep(max(abs(z))/2,2))))
  		text(legend$xtext,legend$ytext,pos=4,
        labels=c(as.character(round(c(0,max(abs(z))/2,max(abs(z)))/factmult,
        legend$digits)),"Positive value","Negative value"))
    }
}

###################       FEEDBACK KRIGING       ####################


.map.shape=function(nodes,region,plots){
    MAP=map('worldHires', region$border,xlim=region$xlim,ylim=region$ylim,col="grey",
        fill=TRUE,plot=plots)
    k=length(MAP)+1
    j1=1
    j=1
    for(i in 1:length(MAP$x)){
    	if(is.na(MAP$x[j])){
            if(j-1-j1>5){
                MAP[[k]]=cbind(MAP$x[j1:(j-1)],MAP$y[j1:(j-1)])
                k=k+1
            }
            j1=j+1
    	}
    	j=j+1
    }
    in.region=0
    for(i in 5:length(MAP)){
    	in.region=in.region+point.in.polygon(nodes[,1],nodes[,2],MAP[[i]][,1],MAP[[i]][,2])
    }
    if(plots){
        points(nodes[in.region>0,],pch=".",col=3)
    }
    list(MAP=MAP,in.region=(in.region>0))
}


krige=function(coordinates,statistic,variog.param,grid,krige.param,plots=TRUE){
    input=list(coordinates=coordinates,statistic=statistic,
        variog.param=variog.param,grid=grid,krige.param=krige.param)
	## projection and scaling of coordinates and statistic
	if(length(grid$proj)>0){
		coord.proj=project(coordinates,proj=grid$proj,degrees=grid$degrees)
		coord.proj=cbind(coord.proj$x,coord.proj$y)
	} else {
		coord.proj=coordinates
	}
	coord.scaling=variog.param$coordinates.scaling
	stat.scaling=variog.param$statistic.scaling
	coord.proj=coord.proj/coord.scaling
	statistic=statistic/stat.scaling
	variog.param$cov.pars=variog.param$cov.pars/c(stat.scaling^2,coord.scaling)
	variog.param$nugget=variog.param$nugget/stat.scaling^2
	## variogram fitting
  	aa=as.geodata(cbind(coord.proj,statistic),coords.col=1:2,data.col=3)
  	maxdist=max(dist(aa$coord))*variog.param$keep.distance
  	cloud1=variog(aa, option = "cloud")
  	bin1=variog(aa,uvec=seq(0,maxdist,l=variog.param$nb.bin),bin.cloud = T)
  	ols.n=variofit(bin1, ini.cov.pars=c(variog.param$cov.pars[1],variog.param$cov.pars[2]),
    	nugget = variog.param$nugget,fix.nugget=variog.param$fix.nugget,
    	weights = "equal")
  	wls.n=variofit(bin1, ini.cov.pars=ols.n$cov.pars,nugget=ols.n$nugget,
    	fix.nugget=variog.param$fix.nugget)
	if(plots){	
		xlab="distance"
		ylab="semivariance"
		if(length(grid$proj)>0){
			xlab=paste(xlab,"in projection space")
		}
		if(coord.scaling!=1){
			xlab=paste("scaled",xlab)
		}
		if(stat.scaling!=1){
			ylab=paste("scaled",ylab)
		}
		plot(bin1,xlab=xlab,ylab=ylab)
  		lines(wls.n, lty = 1, lwd = 2, max.dist = maxdist)
  		## variogram envelopes (under randomness and under the model)
  		env.mc <- variog.mc.env(aa, obj.variog = bin1)
  		env.model <- variog.model.env(aa, obj.variog = bin1,model.pars = wls.n)
  		plot(bin1, envelope = env.mc,xlab=xlab,ylab=ylab)
  		plot(bin1, envelope = env.model,xlab=xlab,ylab=ylab)
  	}
  	## grid shaping
  	gr=expand.grid(grid$x,grid$y)
	gr1=.map.shape(gr,list(border=grid$border,xlim=range(grid$x),ylim=range(grid$y)),
    plots=FALSE)
	MAP=gr1$MAP
	in.region=gr1$in.region
  	## kriging
	if(length(grid$proj)>0){
  		pred.grid=project(gr,proj=grid$proj,degrees=grid$degrees)
  		pred.grid=cbind(pred.grid$x,pred.grid$y)
  	} else {
  		pred.grid=gr
  	}
  	scaled.pred.grid=pred.grid/coord.scaling
  	kc=krige.conv(geodata=aa,locations=scaled.pred.grid,
  		krige=krige.control(obj.model=wls.n,type.krige=krige.param$type.krige,
  			trend.d=krige.param$trend.d))
  	kc$predict=kc$predict*stat.scaling
  	kc$krige.var=kc$krige.var*stat.scaling^2
  	kc$predict[!in.region]=NA
	kc$krige.var[!in.region]=NA
	if(plots){
		boxplot(cbind(kc$predict,sqrt(kc$krige.var)),axes=FALSE,ylab="statistic")
		box()
		axis(1,1:2,labels=c("prediction","sd"))
		axis(2)
		image(kc,loc=gr,col=gray(seq(1,0.1,l=100)),xlab="latitude",ylab="longitude",
        main="kriging prediction")
		contour(kc,loc=gr,add=TRUE)
		map('worldHires', grid$border,xlim=range(grid$x),ylim=range(grid$y),col=1,add=TRUE)	
		image(kc,val=sqrt(kc$krige.var),loc=gr,col=gray(seq(1,0.1,l=100)),
        xlab="latitude",ylab="longitude",main="kriging standard error")
		contour(kc,val=sqrt(kc$krige.var),loc=gr,add=TRUE)
		map('worldHires', grid$border,xlim=range(grid$x),ylim=range(grid$y),col=1,add=TRUE)	
	}
  	list(input=input,in.region=in.region,variofit.wls=wls.n,MAP=MAP,grid=gr,
    krige=kc)
}

krige.test=function(krige.output, subregion, alternative, nb.rand, subregion.coverage=0.8){
	## permute and krige
	coord=krige.output$input$coordinates
	in.region=krige.output$in.region
	grid=krige.output$grid
	in.subregion=(point.in.polygon(grid[,1],grid[,2],subregion$x,subregion$y)>0)
	predict.permute=NULL
	j=0
	k=0
	while(j < nb.rand){
		k=k+1
		translation.index=sample(1:(length(krige.output$krige$predict)-1),1)
		permutation=c((translation.index+1):length(krige.output$krige$predict),
			1:translation.index)
		predict.permute.try=krige.output$krige$predict[permutation]
		if(mean(!is.na(predict.permute.try[in.region & in.subregion]))>subregion.coverage){
			predict.permute=cbind(predict.permute,predict.permute.try)
			j=j+1
			if(j/100==round(j/100)){ 
				print(paste("Number of permutations:",j,"/",nb.rand))
			}
		}
	}
	print(paste("Total number of permutations:",k,"(",k-nb.rand,"permutations led to under-coverage of the subregion under study and were discarded )"))
	## plot
	count.local.nodes=sum(in.region & in.subregion)
	localpred=mean(krige.output$krige$predict[in.region & in.subregion])
	localPRED=as.numeric(colMeans(predict.permute[in.region & in.subregion,],
		na.rm=TRUE))
	if(alternative=="greater"){
		pval=mean(localPRED>localpred)
	} else {
		if(alternative=="less"){
			pval=mean(localPRED<localpred)
		} else {
			stop(paste("Alternative hypothesis",alternative,"is not available."))
		}
	}
    return(new("KT.output",krige.output=krige.output, subregion=subregion,
        averageKrigingPrediction.rand=localPRED,
        averageKrigingPrediction.obs=localpred,
        alternative=alternative, p.value=pval))
}






