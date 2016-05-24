#--------------------------

#' A support function to control the grid. 
#' Primarily used for the GUI in which BcDiag is implemented (RcmdrPlugin.BiclustGUI)

#--------------------------

.checkcurrentgrid <- function(nplots=2,change.x,change.y){
	grid.x <- par()$mfrow[1]
	grid.y <- par()$mfrow[2]
	
	if(grid.x!=1 & grid.y!=1 & grid.x*grid.y>=nplots){}
	else{
		par(mfrow=c(change.x,change.y))
	}
}


#--------------------------

#' A function to transform the isa2 results to biclust results
#' @x, isa object

#--------------------------

isa2biclust <- function(x){
	Parameters <- list(seeddata=x$seeddata,rundata=x$rundata)
	RowxNumber <- (x$rows != 0)
	NumberxCol <- t(x$columns != 0)
	Number <- ncol(x$rows)

	out <- new("Biclust", Parameters=Parameters, RowxNumber=RowxNumber, NumberxCol=NumberxCol, Number=Number)
	return(out)
}


#----------------------

#' a function supports for explore bic function
#' calculates the 'mean','median','variance','quantile' and 'mad'
#' @bic.mat, matrix data either for biclust or outside biclust


#----------------------

exploreCalc<-function(bic.mat){
	#matrix for quantiles
	quantile.25<-matrix(0,nrow(bic.mat),1)
	quantile.50<-matrix(0,nrow(bic.mat),1)
	quantile.75<-matrix(0,nrow(bic.mat),1)
	
	#matrix for mean, median, variance and mad
	mean.res<-matrix(0,nrow(bic.mat),1)
	median.res<-matrix(0,nrow(bic.mat),1)
	var.res<-matrix(0,nrow(bic.mat),1)
	mad.res<-matrix(0,nrow(bic.mat),1)
	
	for(j in 1:nrow(bic.mat)){
		bic.rows<-bic.mat[j,]
		quantile.25[j,]<-quantile(bic.rows,probs=25/100)
		quantile.50[j,]<-quantile(bic.rows,probs=50/100)
		quantile.75[j,]<-quantile(bic.rows,probs=75/100)
		#dimnames(quantile.res)<-list(rownames(bic.mat),names(quantile.res))
		mean.res[j,]<-mean(bic.rows)
		median.res[j,]<-median(bic.rows)
		var.res[j,]<-var(bic.rows)
		mad.res[j,]<-mad(bic.rows)
	}
	dimnames(mean.res)<-list(rownames(bic.mat),"Mean")
	dimnames(median.res)<-list(rownames(bic.mat),"Median")
	dimnames(var.res)<-list(rownames(bic.mat),"Variance")
	dimnames(mad.res)<-list(rownames(bic.mat),"MAD")
	
	dimnames(quantile.25)<-list(rownames(bic.mat),"Bic 25%")
	dimnames(quantile.50)<-list(rownames(bic.mat),"Bic 50%")
	dimnames(quantile.75)<-list(rownames(bic.mat),"Bic 75%")
	quant.res<-list(quantile.25,quantile.50,quantile.75)
	
	#return the calculated summary
	return(list(mean.res,median.res,var.res,mad.res,quant.res))
	
}


#-------------------

#' a function supports for explore bic.
#' used to plot for 'all', 'mean','variance','median', or 'mad'
#' @sbic,@obic; are summary of the bic and outside bicluster respectively
#' @pfor; which plot to figure;'all', 'mean','variance','median', or 'mad'

#-------------------

explorePlot<-function(sbic,obic,pfor=c("all","mean","median","variance","mad","quant")){
	#check for which plot the user is looking for
	pfor<-match.arg(pfor)
	if(pfor=="all"){
		mname<-c("Mean","Median","Variance","MAD")
		#par(mfrow=c(2,2))
		.checkcurrentgrid(4,2,2)
		for(i in 1:4){
			plot(c(sbic[[i]],obic[[i]]),type="n",col=4,ylab="",xlab="",main=mname[i])
			lines(c(sbic[[i]],obic[[i]]),type="l",col=1)
			lines(sbic[[i]],col=2)
		}
	}
	#par(mfrow=c(1,1))
	if(pfor=="quant"){
		#par(mfrow=c(1,1))
		bquant1<-sbic[[5]][[1]];oquant1<-obic[[5]][[1]]
		bquant2<-sbic[[5]][[2]];oquant2<-obic[[5]][[2]]
		bquant3<-sbic[[5]][[3]];oquant3<-obic[[5]][[3]]
		
		plot(c(bquant3,oquant3),lty=1,type="n",lwd=1,xlab="",ylab="")
		lines(c(bquant3,oquant3),lty=1,type="l",col=2,lwd=1)
		lines(c(bquant2,oquant2),lty=1,type="l",col=3,lwd=1)
		lines(c(bquant1,oquant1),lty=1,type="l",col=4,lwd=1)
		legend("topleft",c("75%","50%","25%"),col=c(2,3,4),lty=c(1,1,1))
	}
	if(pfor=="mean"){  
		plot(c(sbic[[1]],obic[[1]]),type="n",xlab="",ylab="",main=pfor)
		lines(c(sbic[[1]],obic[[1]]),type="l",lty=1,lwd=1)
		lines(sbic[[1]],col=2,lty=1)
	}
	if(pfor=="median"){
		plot(c(sbic[[2]],obic[[2]]),type="n",xlab="",ylab="",main=pfor)
		lines(c(sbic[[2]],obic[[2]]),type="l",lty=1,lwd=1)
		lines(sbic[[2]],col=2,lty=1)
	}
	if(pfor=="variance"){
		plot(c(sbic[[3]],obic[[3]]),type="n",xlab="",ylab="",main=pfor)
		lines(c(sbic[[3]],obic[[3]]),type="l",lty=1,lwd=1)
		lines(sbic[[3]],col=2,lty=1)
	}
	if(pfor=="mad"){
		plot(c(sbic[[4]],obic[[4]]),type="n",xlab="",ylab="",main=pfor)
		lines(c(sbic[[4]],obic[[4]]),type="l",lty=1,lwd=1)
		lines(sbic[[4]],col=2,lty=1)
	}
}


#----------------------

#' a supportive function for profileplot.bic,
#' indexedBic: a function to check the name of the method
#' and returns the list parameter of the required indcies of the biclust
#' @ bres, the biclust object
#' @ mname, name of the method to be applied for the biclust
#' @ dset,bres,mname,bnum has similar explanation as summary.bic function
#' outcome : returns the two indcies; the indg and indc
#' @ indg; index for the biclust genes.
#' @ indc; index for the biclust conditions.

#----------------------
indexedBic<-function(dset,bres,mname=c("fabia","isa2","biclust","bicare"),bnum,fabia.thresZ=0.5,fabia.thresL=NULL){
	# which biclust object is it; 
	check<-match.arg(mname)
	l<-bnum
	if(check=="fabia"){
		#Extract biclusters:
		#get the biclust index inside the dset 
		resf <- extractBic(bres,thresZ=fabia.thresZ,thresL=fabia.thresL)
		bg<-resf$numn[l,]$numng
		bc<-resf$numn[l,]$numnp
		# the two indecies	
		indg<-bg
		indc<-bc
	}
	if(check=="isa2"){
		#convert to biclust and get the biclust indecies
		resi<-isa2biclust(bres)
		indg<-which(resi@RowxNumber[,l])
		indc<-which(resi@NumberxCol[l,])
		
	}
	if(check=="biclust"){
		indg<-which(bres@RowxNumber[,l])
		indc<-which(bres@NumberxCol[l,])
		
	}
	if(check=="bicare"){
		x <- bres
		Parameters <- list(numberofbicluster=x$param[1,2],residuthreshold=x$param[2,2],genesinitialprobability=x$param[3,2],samplesinitialprobability=x$param[4,2],numberofiterations=x$param[5,2],date=x$param[6,2])
		RowxNumber <- t(x$bicRow==1)
		NumberxCol <- x$bicCol==1
		Number <- as.numeric(dim(RowxNumber)[2])
		info <- list()
		resbi <- new("Biclust",Parameters=Parameters,RowxNumber=RowxNumber,NumberxCol=NumberxCol,Number=Number,info=info)
		indg<-which(resbi@RowxNumber[,l])
		indc<-which(resbi@NumberxCol[l,])
	}
	
	# Warning for bicluster with 1 row and 1 column
	if(length(indg)==1 & length(indc)==1){stop("Bicluster only has 1 row and 1 column",call.=FALSE)}

	return(list(indg,indc))
}

#----------------

#' A function to plot summaries only for biclust data.
#' @dset stands for the dataset,
#' @bres stands for biclust result, should be one of fabia, biclust or isa2
#' @fit a parameter takes string; 'all'(all plots in one frame), 'mean', 'median', 'variance' or 'mad'.
#' @gby, group by 'conditions' or 'genes'

#-----------------

plotOnlybic<-function(ball,fit="all",gby){
	
	if(fit=="all"){
		#biclust genes
		sname<-c("Median","Mean","Variance","MAD")
		#par(mfrow=c(2,2))
		.checkcurrentgrid(4,2,2)
		for(i in 1:length(ball)){
			plot(ball[[i]],type="n",main=sname[i],ylab="",xlab=gby)
			lines(ball[[i]],type="l",col=i+1)
			
		}
	}
	#par(mfrow=c(1,1))
	if(fit=="median")
	{
		plot(ball[[1]],type="n",main="Median",ylab="",xlab="Condtions")
		lines(ball[[1]],type="l",col=2)
	}
	if(fit=="mean")
	{
		plot(ball[[2]],type="n",main="Mean",ylab="",xlab=gby)
		lines(ball[[2]],type="l",col=3)
		
	}
	if(fit=="variance")
	{
		plot(ball[[3]],type="n",main="Variance",ylab="",xlab=gby)
		lines(ball[[3]],type="l",col=4)
		
	}
	if(fit=="mad")
	{
		plot(ball[[4]],type="n",main="MAD",ylab="",xlab=gby)
		lines(ball[[4]],type="l",col=5)
		
	}
}


#-----------------

#' supportive function for profileBic
#' plots profile in four different types of plots.
#' outcome: four plots in one: profile lines,boxplot,histogram and 3d plots.
#' @ indg; index for the biclust genes.
#' @ indc; index for the biclust conditions.
#' @gby; group by 'genes' or 'conditions'

#------------------
profileAll<-function(dset,indg,indc,grp,gby="genes",teta=120,ph=30){
	
	if(gby=="conditions"){
		
		#group the genes in to two.
		cnams <- colnames(dset)
		grp <- rep(1, length(cnams))
		grp[indc] <- 2
		d<-dset[indg, order(grp, decreasing=TRUE)]
		dbc<-dset[indg,(grp==2)]		
		#lines 
		matplot(y =t(d),type ="n",col="green3", xlab="Condtions",ylab="Expression", axes=T, pch = rep(1, ncol(d)))
		matlines(y = t(d), type = "l",lty = rep(1, nrow(d)) ,col="green3", lwd = 1, pch = 1)
		matlines(y = t(dbc), type = "l",lty = rep(1, nrow(d)) ,col=2, lwd = 1, pch = 1)
		axis(2)
		box()
		
		#boxplot 
		boxplot.matrix(d,col="green3",main="",xlab="",axes=T,lty=1)
		boxplot.matrix(dbc,add=T,col="red",axes=F,xlab="",lty=1)
		box()
		
		#hist 
		hist(d,col="green3",xlab="",main="",lty=1)
		hist( dbc,col="red",add=T,xlab="",main="",lty=1)
		box()
		
		#3D 
		d1<-c(1:nrow(d))
		d2<-c(1:ncol(d))
		
		fill <- matrix("green3", nrow = nrow(d), ncol = ncol(d))
		fill[,(grp==2)] <- "red";fill<-sort(fill, decreasing=T)
		persp(d1,d2,d,theta = teta, phi = ph, expand = 0.5, col = fill,
				ltheta = 120, shade = 0.75, ticktype = "detailed",
				,xlab="Genes",ylab="Condtions",zlab="Gene expression")
		box()
	}
	if(gby=="genes"){
		
		#group the genes in to two.
		rnams <- rownames(dset)
		grp <- rep(1, length(rnams))
		grp[indg] <- 2
		d<-dset[order(grp, decreasing=TRUE),indc]	
		dbc<-dset[(grp==2),indc]
		
		#lines
		matplot(y = d, type = "n",xlab="",ylab="Expression", axes=T, pch = rep(1, ncol(dset)))
		matlines(y = d, type = "l", col="black",lty = rep(1, nrow(dset)),lwd = 1, pch = 1) 
		matlines(y = dbc, type = "l",col="red",lty = rep(1, nrow(dset)) , lwd = 1, pch = 1)
		box()
		
		#boxplot
		boxplot.matrix(t(d),col="black",main="",xlab="",axes=T)
		lines(t(dbc),type="l",col="red")
		box()
		
		#histogram
		hist(d,col="black",xlab="",main="")
		hist(dbc,col="red",add=T,xlab="",main="")
		box()
		
		#3D
		d1<-c(1:nrow(d))
		d2<-c(1:ncol(d))
		fill <- matrix("black", nrow = nrow(d), ncol = ncol(d))
		fill[(grp==2),] <- "red";fill<-sort(fill,decreasing=T)
		persp(d1,d2,d,theta = teta, phi = ph, expand = 0.5, col = fill,
				ltheta = 120, shade = 0.75, ticktype = "detailed",
				,xlab="Genes",ylab="Condtions",zlab="Gene expression")
		box()
		
	}
	
}



