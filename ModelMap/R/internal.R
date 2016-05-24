
#############################################################################################
#############################################################################################
########################## Initialize graphics device #######################################
#############################################################################################
#############################################################################################


initialize.device<-function(	PLOTfn,
					DEVICE.TYPE,
					res,
					device.width,
					device.height,
					units,
					pointsize,
					cex
){
	if(DEVICE.TYPE=="default"){dev.new(width = device.width, height = device.height, pointsize = pointsize)}
	if(DEVICE.TYPE=="jpeg"){jpeg(filename=paste(PLOTfn,".jpg",sep=""),width = device.width, height = device.height, res=res, units=units, pointsize=pointsize)}
	if(DEVICE.TYPE=="pdf"){pdf(file=paste(PLOTfn,".pdf",sep=""),width = device.width, height = device.height, pointsize=pointsize)}
	if(DEVICE.TYPE=="png"){png(filename=paste(PLOTfn,".png",sep=""),width = device.width, height = device.height, res=res, units=units, pointsize=pointsize)}
	if(DEVICE.TYPE=="postscript"){postscript(file=paste(PLOTfn,".ps",sep=""),width = device.width, height = device.height, pointsize=pointsize)}
	if(DEVICE.TYPE=="tiff"){tiff(filename=paste(PLOTfn,".tif",sep=""),width = device.width, height = device.height, res=res, units=units, pointsize=pointsize)}	
}




#############################################################################################
#############################################################################################
###################################### Correlation ##########################################
#############################################################################################
#############################################################################################


correlation.function<-function(	qdata,
						predList,
						predFactor,
						MODELpredfn,
						main=basename(MODELpredfn),

						device.type="default",
						res=72,
						device.width=7,
						device.height=7,
						units="in",
						pointsize=12,
						cex=par()$cex
						){

if(!"none"%in%device.type){
	for(i in 1:length(device.type)){
		#print(paste("Device Type:",device.type[i]))

		### Output filenames ###

		CORRPLOTfn<-paste(MODELpredfn,"_corrplot",sep="")

		initialize.device(	PLOTfn=CORRPLOTfn,DEVICE.TYPE=device.type[i],
						res=res,device.width=device.width,device.height=device.height,
						units=units,pointsize=pointsize,cex=cex)

		opar<-par(cex=cex,oma=c(0,0,3,0))
		
		M<-cor(qdata[,predList[!predList%in%predFactor]])
		corrplot(M,order="AOE",type="upper",tl.pos="tp")
		corrplot(M,add=TRUE, type="lower", method="number",order="AOE", col="black",
				diag=FALSE,tl.pos="n", cl.pos="n")
		#corrplot(M)
		mtext(main,side=3,line=1,cex=1.3*cex,outer=TRUE)
		mtext("Predictor Correlation",side=3,line=-0.5,cex=1.3*cex,outer=TRUE)

		par(opar)
		if(!device.type[i]%in%c("default","none")){dev.off()}
	}
}
}



#############################################################################################
#############################################################################################
########################## model.explore - Mask Function ####################################
#############################################################################################
#############################################################################################

mask.fun <- function(RB, PRED.RANGE, predContinuous, predFactor, N.Continuous, N.Factor, filename) {

	print(paste("nlayers(RB) =",nlayers(RB)))

	MB <- brick(RB[[1]],nl=nlayers(RB),values=FALSE, filename=filename, overwrite=TRUE)
	names(MB)<-names(RB)
 	filename <- trim(filename)

	print("created MB")
	
# 	if (canProcessInMemory(MB, 3)) {
#		v <- getValues(RB)
#		
#		if(N.Continuous>0){for(pred in predContinuous){
#			v[,pred][v[,pred]<PRED.RANGE[[pred]][1]]<-NA
#			v[,pred][v[,pred]>PRED.RANGE[[pred]][2]]<-NA
#		}}
#		if(N.Factor>0){for(pred in predFactor){
#			v[,pred][!v[,pred]%in%PRED.RANGE[[pred]]]<-NA
#		}}
#		MB<- setValues(MB, v, filename=filename, overwrite=TRUE)
#	}else{
		bs <- blockSize(MB)
		MB <- writeStart(MB, filename, overwrite=TRUE)
		for (i in 1:bs$n) {
				v <- as.matrix(getValues(RB, row=bs$row[i], nrows=bs$nrows[i] ))
				if(nlayers(RB)==1){colnames(v)<-names(PRED.RANGE)}
				if(N.Continuous>0){for(pred in predContinuous){
					v[,pred][v[,pred]<PRED.RANGE[[pred]][1]]<-NA
					v[,pred][v[,pred]>PRED.RANGE[[pred]][2]]<-NA
				}}
				if(N.Factor>0){for(pred in predFactor){
					v[,pred][!v[,pred]%in%PRED.RANGE[[pred]]]<-NA
				}}
				MB <- writeValues(MB, v, bs$row[i])	
		}
		MB <- writeStop(MB)
#	}
 	return(MB)
 }

##############################################

mask.fun <- function(RB, PRED.RANGE, predContinuous, predFactor, N.Continuous, N.Factor, filename) {

	print(paste("nlayers(RB) =",nlayers(RB)))

	MB <- brick(RB[[1]],nl=nlayers(RB),values=FALSE, filename=filename, overwrite=TRUE)
	names(MB)<-names(RB)
 	filename <- trim(filename)

	print("created MB")
	
# 	if (canProcessInMemory(MB, 3)) {
#		v <- getValues(RB)
#		
#		if(N.Continuous>0){for(pred in predContinuous){
#			v[,pred][v[,pred]<PRED.RANGE[[pred]][1]]<-NA
#			v[,pred][v[,pred]>PRED.RANGE[[pred]][2]]<-NA
#		}}
#		if(N.Factor>0){for(pred in predFactor){
#			v[,pred][!v[,pred]%in%PRED.RANGE[[pred]]]<-NA
#		}}
#		MB<- setValues(MB, v, filename=filename, overwrite=TRUE)
#	}else{
		bs <- blockSize(MB)
		MB <- writeStart(MB, filename, overwrite=TRUE)
		for (i in 1:bs$n) {
				v <- as.matrix(getValues(RB, row=bs$row[i], nrows=bs$nrows[i] ))
				if(nlayers(RB)==1){colnames(v)<-names(PRED.RANGE)}
				if(N.Continuous>0){for(pred in predContinuous){
					v[,pred][v[,pred]<PRED.RANGE[[pred]][1]]<-NA
					v[,pred][v[,pred]>PRED.RANGE[[pred]][2]]<-NA
				}}
				if(N.Factor>0){for(pred in predFactor){
					v[,pred][!v[,pred]%in%PRED.RANGE[[pred]]]<-NA
				}}
				MB <- writeValues(MB, v, bs$row[i])	
		}
		MB <- writeStop(MB)
#	}
 	return(MB)
 }


#############################################################################################
#############################################################################################
########################## model.explore - Mask Graphic ####################################
#############################################################################################
#############################################################################################

mask.graphic<-function(	RL=RL,
				ML=ML,
				#MAXCELL=100000,
				OUTfn,
				main="Extrapolation Mask - All Predictors",
				device.type="default",
				res=72,
				device.width=0.9*device.height,
				device.height=7,
				units="in",
				pointsize=12,
				cex=par$cex){

#if(ncell(RL)>MAXCELL){RL<-sampleRegular(RL, size=MAXCELL, asRaster=TRUE)}
#if(ncell(ML)>MAXCELL){ML<-sampleRegular(ML, size=MAXCELL, asRaster=TRUE)}

if(!"none"%in%device.type){
	for(i in 1:length(device.type)){

		print(paste("device.type =",device.type[i]))

		initialize.device(	PLOTfn=OUTfn,DEVICE.TYPE=device.type[i],
						res=res,device.width=device.width,device.height=device.height,
						units=units,pointsize=pointsize,cex=cex)

		###img.rast: 0=NA(grey) 1=MASK(black) 2=DATA(green)
		img.rast<-RL+1
		img.rast[is.na(ML)]<-1
		img.rast[is.na(RL)]<-0
		
		opar<-par(cex=cex,mar=c(4,3,4,1)+.1)
		image( img.rast, 
		       col = c("grey90","black","palegreen2"),
			 zlim=c(0,2),
		       xlab="", ylab="", 
		       asp=1, bty="n", main="")

#		if(cellStats(is.na(RL),sum)>0){
#			NAval.RL<-is.na(RL)
#			NAval.RL[!NAval.RL]<-NA
#			image( NAval.RL, 
#			       col = "grey90",
#			       xlab="", ylab="", 
#			       #zlim=zlim,
#			       #asp=1, bty="n", main="",
#				 add=TRUE)
#		}

		legend( #x=xmax(mapgrid),y=ymax(mapgrid),
			  x="topright",
			  bg="white",
		        legend=c("Data","NA","MASK"),
		        fill=c("palegreen2","grey90","black"),
		        #bty="n",
			  #title=pred,
		        cex=cex)
		mtext(main,side=3,line=2.2,cex=1.5*cex)

		par(opar)

		if(!device.type[i]%in%c("default","none")){dev.off()}

	}
}
}



#############################################################################################
#############################################################################################
################## model.explore - Categorical Predictors - Three plots #####################
#############################################################################################
#############################################################################################


explore.categorical<-function(qdata,
					response.name,
					response.type,
					response.colors,
					pred,
					pred.range,
					rast.range,
					OUTfn,
					main=basename(OUTfn),

					pred.rast,
					pred.mask,
					#MAXCELL=100000,
	
					col.cat,

					device.type="default",
					res=72,
					device.width=device.height*1.8,
					device.height=7,
					units="in",
					pointsize=12,
					cex=par$cex
					){

		LEV.train  <-pred.range
		LEV.rast   <-rast.range
		LEV.all    <-sort(unique(c(LEV.train,LEV.rast)))
		LEV.missing<-LEV.rast[!LEV.rast%in%LEV.train]

		#cex.names <- if(length(LEV.all)<=10){0.8*cex}else{0.6*cex}
		cex.names  <- if(length(LEV.all)<=10){0.8*cex}else{if(length(LEV.all)<=25){0.6*cex}else{0.4*cex}}
		cex.legend <- if(length(LEV.all)<=10){1.0*cex}else{if(length(LEV.all)<=25){0.8*cex}else{0.6*cex}}

		legend.label<-LEV.all
		#if(any(!LEV.rast%in%LEV.train)){
		#	legend.label[!LEV.rast%in%LEV.train]<-paste("**",legend.label[!LEV.rast%in%LEV.train],"**",sep="")
		#}
		if(any(!LEV.rast%in%LEV.train)){
			legend.label[legend.label%in%LEV.missing]<-paste("**",legend.label[legend.label%in%LEV.missing],"**",sep="")
		}
		Nlab<-length(legend.label)
		legend.colors<-rep(col.cat,length.out=Nlab)

		### ### ### SCREEN 4 ### ### ###
		#set screen to left half

		###moved sampling to model.explore()###
		#pred.stack<-stack(list(pred=pred.rast,mask=pred.mask))
		#
		#if(ncell(pred.rast)>MAXCELL){
		#	samp.stack<-sampleRegular(pred.stack, size=MAXCELL, asRaster=TRUE)
		#}else{samp.stack<-pred.stack}
		#
		#samp.rast<-samp.stack[["pred"]]
		#samp.mask<-mask(samp.stack[["pred"]],samp.stack[["mask"]])

		#img.data = data raster with pixels outside of range of training masked to NA#
		img.data<-mask(pred.rast,pred.mask)

		#imp.mask = 0 if original raster is NAval, and 1 if original raster is not NAval, irregardless of training range#
		img.mask<-pred.rast
		img.mask[is.na(pred.rast)]<-0
		img.mask[!is.na(pred.rast)]<-1

		rast.lab<-unique(img.data)
		Nrast<-length(rast.lab)
		v <- getValues(img.data)
		v <- (1:Nrast)[match(v,unique(img.data))]
		img.data <- setValues(pred.mask, v)

		### ### ### SCREEN 1 ### ### ###
		#set screen to right upper
		if(response.type=="continuous"){
			means.resp<-tapply(qdata[,response.name],qdata[,pred],mean)}

		if(response.type=="binary"){
			counts.resp<-table(qdata[,response.name],qdata[,pred])
			counts.resp<-counts.resp[c(2,1),]
			counts.resp.stacked<-apply(counts.resp,2,sum)
			counts.resp.relative<-counts.resp/matrix(counts.resp.stacked,nrow(counts.resp),ncol(counts.resp),byrow=TRUE)}


		if(response.type=="categorical"){
			counts.resp<-table(qdata[,response.name],qdata[,pred])
			counts.resp.stacked<-apply(counts.resp,2,sum)
			counts.resp.relative<-counts.resp/matrix(counts.resp.stacked,nrow(counts.resp),ncol(counts.resp),byrow=TRUE)}
		
		### ### ### SCREEN 2 ### ### ###
		#set screen to right middle
		TABLE.train<-table(qdata[,pred])
		counts.train<-rep(0,Nlab)
		names(counts.train)<-LEV.all
		counts.train[names(TABLE.train)]<-TABLE.train

		### ### ### SCREEN 3 ### ### ###
		#set screen to right lower
		TABLE.rast<-freq(pred.rast)
		counts.rast<-rep(0,Nlab)
		names(counts.rast)<-LEV.all
		counts.rast[as.character(TABLE.rast[,1][!is.na(TABLE.rast[,1])])]<-TABLE.rast[,2][!is.na(TABLE.rast[,1])]


if(!"none"%in%device.type){
	for(i in 1:length(device.type)){
		#print(paste("Device Type:",device.type[i]))

		### Output filenames ###

		initialize.device(	PLOTfn=OUTfn,DEVICE.TYPE=device.type[i],
						res=res,device.width=device.width,device.height=device.height,
						units=units,pointsize=pointsize,cex=cex)

		layout(matrix(c(4,1,4,2,4,3),3,2,byrow=TRUE))

		### ### ### SCREEN 1 ### ### ###
		#set screen to right upper

		if(response.type=="continuous"){
			opar <- par(cex=cex, mgp=c(1.5, .5, 0),mar=c(3,3,2,1)+.1)
			barplot(	means.resp, 
					xlab=pred, ylab=response.name,
					main="",
					cex.names=cex.names,
					col=legend.colors[legend.label%in%names(means.resp)])
			mtext("Mean Response - Training Data", side=3, line=0.5, cex=1.2*cex, font=2 )
			par(opar)}

		if(response.type=="binary"){
			opar <- par(cex=cex, mgp=c(1.5, .5, 0),mar=c(3,3,2,5)+.1, xpd=TRUE)
			color=c("gray80","gray20")
			#color=c("red","blue")
			RelBP<-TRUE
			if(!RelBP){
			BP<-barplot(counts.resp,
					col=color,
					names.arg=colnames(counts.resp),cex.names=cex.names,
					xlab=pred, ylab="frequency")
			}else{
			BP<-barplot(counts.resp.relative,
					col=color,
					names.arg=colnames(counts.resp),cex.names=cex.names,
					xlab=pred, ylab="proportion")
			}
			mtext("Training Data", side=3, line=0.5, cex=1.2*cex, font=2 )
			#axis(1,at=seq(min(BP),max(BP),5))
			legend(	x="topright",inset=c(-.2,-.2),bg="white",title=response.name,
					legend=c("absent","present"),fill=rev(color),cex=0.8*cex)
			par(opar)
		}

		if(response.type=="categorical"){
			opar <- par(cex=cex, mgp=c(1.5, .5, 0),mar=c(3,3,2,5)+.1, xpd=TRUE)
			RelBP<-TRUE
		   if(!RelBP){
			BP<-barplot(counts.resp,
					col=response.colors$colors,
					names.arg=colnames(counts.resp),cex.names=cex.names,
					xlab=pred, ylab="frequency")
			mtext("Response - Training Data", side=3, line=0.5, cex=1.2*cex, font=2 )
			#axis(1,at=seq(min(BP),max(BP),5))
			legend(	x="topright",inset=c(-.2,-.2),title=response.name,bg="white",
					legend=rev(response.colors$category),fill=rev(response.colors$colors),cex=0.8*cex)}

		   if(RelBP){
			BP<-barplot(counts.resp.relative,
					col=response.colors$colors,
					names.arg=colnames(counts.resp.relative),cex.names=cex.names,
					xlab=pred, ylab="proportion")
			mtext("Response - Training Data", side=3, line=0.5, cex=1.2*cex, font=2 )
			#axis(1,at=seq(min(BP),max(BP),5))
			legend(	x="topright",inset=c(-.2,-.2),title=response.name,bg="white",
					legend=rev(response.colors$category),fill=rev(response.colors$colors),cex=0.8*cex)}

			par(opar)
		}
		
		
		### ### ### SCREEN 2 ### ### ###
		#set screen to right middle
		opar <- par(cex=cex, mgp=c(1.5, .5, 0),mar=c(3,3,2,1)+.1)
		barplot(counts.train, xlab=pred, ylab="frequency", main="",cex.names=cex.names,col=legend.colors)
		mtext("Number of Plots - Training", side=3, line=0.5, cex=1.2*cex, font=2 )
		par(opar)

		### ### ### SCREEN 3 ### ### ###
		#set screen to right lower 
		opar <- par(cex=cex, mgp=c(1.5, .5, 0),mar=c(3,3,2,1)+.1)
		barplot(counts.rast, xlab=pred, ylab="frequency", main="",cex.names=cex.names,col=legend.colors)
		mtext("Number of Pixels - Raster", side=3, line=0.5, cex=1.2*cex, font=2 )

		par(opar)

		### ### ### SCREEN 4 ### ### ###
		#set screen to left half
		
		opar<-par(cex=cex,mar=c(4,3,4,1)+.1)

		image( img.mask, 
		       col = c("grey80","black"),
		       xlab="", ylab="", 
		       zlim=c(0,1),
		       asp=1, bty="n", main="")

		image( img.data, 
		       col = legend.colors[match(rast.lab,legend.label)],
		       #xlab="", ylab="", 
		       #asp=1, bty="n", main="",
			 add=TRUE)

		legend( #x=xmax(mapgrid),y=ymax(mapgrid),
			  x="topright",
			  bg="white",
		        legend=c(legend.label,"NA","MASK"),
		        fill=c(legend.colors,"grey90","black"),
		        #bty="n",
			  #title=pred,
		        cex=cex.legend)
		mtext(paste("Extrapolation Mask -",pred),side=3,line=2.2,cex=1.5*cex)
		mtext("(masked to categories of training)",line=1,cex=cex)
		if(any(!LEV.rast%in%LEV.train)){
			mtext("** = predictor level not found in training data, thus masked out of map",
				side=1, line=3,cex=cex,adj=0)}

		par(opar)

		if(!device.type[i]%in%c("default","none")){dev.off()}
	}
}
}


#############################################################################################
#############################################################################################
######################## model.explore - Continuous Predictors ##############################
#############################################################################################
#############################################################################################


explore.continuous<-function(qdata,
					response.name,
					response.type,
					response.colors,
					pred,
					pred.range,
					OUTfn,
					main=basename(OUTfn),

					pred.rast,
					pred.mask,
					#MAXCELL=100000,

					col.ramp,
					device.type="default",
					res=72,
					device.width=device.height*1.8,
					device.height=7,
					units="in",
					pointsize=12,
					cex=par$cex
					){

		lim<-range(range(qdata[,pred]),minValue(pred.rast),maxValue(pred.rast))

		###moved sampling to model.explore()###
		#pred.stack<-stack(list(pred=pred.rast,mask=pred.mask))
		#
		#if(ncell(pred.rast)>MAXCELL){
		#	samp.stack<-sampleRegular(pred.stack, size=MAXCELL, asRaster=TRUE)
		#}else{samp.stack<-pred.stack}
		#
		#samp.rast<-samp.stack[["pred"]]
		#samp.mask<-mask(samp.stack[["pred"]],samp.stack[["mask"]])

		#img.data = data raster with pixels outside of range of training masked to NA#
		img.data<-mask(pred.rast,pred.mask)

		#imp.mask = 0 if original raster is NAval, and 1 if original raster is not NAval, irregardless of training range#
		img.mask<-pred.rast
		img.mask[is.na(pred.rast)]<-0
		img.mask[!is.na(pred.rast)]<-1

		HIST.qdata<-hist(qdata[,pred], plot=FALSE)
		HIST.rast<-hist(pred.rast,plot=FALSE)

		### ### ### SCREEN 5 ### ### ###
		#set screen to left half
		#print("Screen 5 calc")

		TRAIN.min<-min(qdata[,pred])
		TRAIN.max<-max(qdata[,pred])
		RAST.min <-minValue(pred.rast)
		RAST.max <-maxValue(pred.rast)

		zlim<-range(TRAIN.min,TRAIN.max,RAST.min,RAST.max)
		legend.label<-rev(pretty(zlim,n=5))
		Ncol<-length(col.ramp)
		Nlab<-length(legend.label)
		legend.colors<-col.ramp[rev(round(  ((((1:Nlab)-1)/(Nlab-1))*(Ncol-1))+1  ))]

		### ### ### SCREEN 1 ### ### ###
		#set screen to right half. upper left
		#print("Screen 1 calc")

		if(response.type=="continuous"){
			#lm.pred<-lm(as.formula(paste(response.name,"~",pred)),data=qdata)
			gf<-as.formula(paste(response.name,"~s(",pred,")",sep=""))
			gam.pred<-mgcv::gam(gf,family=gaussian(link=identity),data=qdata)
			gam.line<-data.frame(seq(min(qdata[,pred],na.rm=TRUE),max(qdata[,pred],na.rm=TRUE),length = 100))
			names(gam.line)<-pred
			gam.line[,response.name]<-mgcv::predict.gam(gam.pred, newdata = gam.line)
			ylim=range(qdata[,response.name])
			COR.pearson<-round(cor(qdata[,pred],qdata[,response.name],method="pearson"),2)
			#COR.spearman<-round(cor(qdata[,pred],qdata[,response.name],method="spearman"),2)
		}

		if(response.type=="binary"){
			
			counts.resp<-rbind(	hist(qdata[qdata[,response.name]==1,pred],breaks=HIST.qdata$breaks,plot=FALSE)$counts,
							hist(qdata[qdata[,response.name]==0,pred],breaks=HIST.qdata$breaks,plot=FALSE)$counts)
			counts.resp.stacked<-apply(counts.resp,2,sum)
			counts.resp.relative<-counts.resp/matrix(counts.resp.stacked,nrow(counts.resp),ncol(counts.resp),byrow=TRUE)
			}

		if(response.type=="categorical"){

			Ncat<-nrow(response.colors)
			counts.resp<-hist(qdata[qdata[,response.name]==response.colors$category[1],pred],breaks=HIST.qdata$breaks,plot=FALSE)$counts
			if(Ncat>1){for(j in 2:Ncat){
				counts.resp<-rbind(	counts.resp,
								hist(qdata[qdata[,response.name]==response.colors$category[j],pred],breaks=HIST.qdata$breaks,plot=FALSE)$counts)
			}}
			counts.resp.stacked<-apply(counts.resp,2,sum)
			counts.resp.relative<-counts.resp/matrix(counts.resp.stacked,nrow(counts.resp),ncol(counts.resp),byrow=TRUE)

			#counts.resp<-tapply(qdata[,pred],qdata[,response.name],function(x,breaks){hist(x,breaks=breaks,plot=FALSE)$counts},simplify=TRUE,breaks=HIST.qdata$breaks)
			}

		### ### ### SCREEN 2 ### ### ###
		#set screen to right half. upper right
		#print("Screen 2 calc")

		bp.data<-boxplot(qdata[,pred],plot=FALSE)
		bp.rast<-boxplot(pred.rast,plot=FALSE)


		bp<-bp.data
		bp$stats <- cbind(bp$stats,bp.rast$stats)
		bp$n     <- c(bp$n,bp.rast$n)
		bp$conf  <- cbind(bp$conf,bp.rast$conf)
		bp$names <- c("Train","Raster")	

		### ### ### SCREEN 3 ### ### ###
		#set screen to right half. lower left
		#print("Screen 3 calc")

		MIDS<-HIST.qdata$mids
		MIDS[MIDS<zlim[1]]<-zlim[1]
		MIDS[MIDS>zlim[2]]<-zlim[2]

		#print("zlim")
		#print(zlim)
		#print("qdata mids")
		#print(MIDS)

		hcol.qdata<-col.ramp[round(((MIDS-zlim[1])/(zlim[2]-zlim[1]))*(Ncol-1))+1]	

		### ### ### SCREEN 4 ### ### ###
		#set screen to right half. lower right
		#print("Screen 4 calc")

		MIDS<-HIST.rast$mids
		MIDS[MIDS<zlim[1]]<-zlim[1]
		MIDS[MIDS>zlim[2]]<-zlim[2]

		#print("zlim")
		#print(zlim)
		#print("rast mids")
		#print(MIDS)

		hcol.rast<-col.ramp[round(((MIDS-zlim[1])/(zlim[2]-zlim[1]))*(Ncol-1))+1]

if(!"none"%in%device.type){
	for(i in 1:length(device.type)){
		#print(paste("Device Type:",device.type[i]))

		### Output filenames ###

		initialize.device(	PLOTfn=OUTfn,DEVICE.TYPE=device.type[i],
						res=res,device.width=device.width,device.height=device.height,
						units=units,pointsize=pointsize,cex=cex)

		#layout(matrix(c(1,2,5,5,3,4,5,5),2,4,byrow=TRUE))
		layout(matrix(c(5,5,1,2,5,5,3,4),2,4,byrow=TRUE))

		
		### ### ### SCREEN 1 ### ### ###
		#set screen to right half. upper left
		#print("Screen 1 print")

		if(response.type=="continuous"){
			#par(new=TRUE)
			#plot(gam.pred,ylim=range(qdata[,response.name]),rug=FALSE,col="red")

			opar <- par(cex=cex, mgp=c(1.5, .5, 0),mar=c(3,3,4,1)+.1)
			plot(qdata[,pred],qdata[,response.name],xlab=pred,ylab=response.name,main="",ylim=ylim)
			#abline(lm.pred)
			lines(gam.line[,pred],gam.line[,response.name],col="red",lwd=2)

			mtext("Training Data", side=3, line=2.5, cex=1.25*cex, font=2 )
			mtext(paste(response.name,"~s(",pred,")",sep=""),side=3,line=0.3,cex=cex)
			mtext(paste("  corr ",COR.pearson," ",sep=""),side=3,line=-1,adj=0,cex=0.8*cex)
			par(opar)
		}

		if(response.type=="binary"){
			opar <- par(cex=cex, mgp=c(1.5, .5, 0),mar=c(3,3,5,3)+.1,xpd=NA)
			color=c("gray80","gray20")
			#color=c("red","blue")
			RelBP<-TRUE
			if(!RelBP){
			  BP<-barplot(counts.resp,
					col=color,
					names.arg=HIST.qdata$mids,cex.names=0.8*cex,
					xlab=pred, ylab="frequency")
			}else{
			  BP<-barplot(counts.resp.relative,
					col=color,
					names.arg=HIST.qdata$mids,cex.names=0.8*cex,
					xlab=pred, ylab="proportion")
			}
			mtext("Training Data", side=3, line=2.5, cex=1.25*cex, font=2 )
			#axis(1,at=seq(min(BP),max(BP),5))
			legend(	x="topright",inset=c(-.2,-.2),bg="white",title=response.name,
					legend=c("absent","present"),fill=rev(color),cex=0.8*cex)
			par(opar)
		}

		if(response.type=="categorical"){
			opar <- par(cex=cex, mgp=c(1.5, .5, 0),mar=c(3,3,5,3)+.1,xpd=NA)

			RelBP<-TRUE
			if(!RelBP){
			BP<-barplot(counts.resp,
					col=response.colors$colors,
					names.arg=HIST.qdata$mids,cex.names=0.8*cex,
					xlab=pred, ylab="frequency")
			mtext("Training Data", side=3, line=2.5, cex=1.25*cex, font=2 )
			#axis(1,at=seq(min(BP),max(BP),5))
			legend(	x="topright",inset=c(-.3,-.2),title=response.name,bg="white",
					legend=rev(response.colors$category),fill=rev(response.colors$colors),cex=0.6*cex)}

			if(RelBP){
			BP<-barplot(counts.resp.relative,
					col=response.colors$colors,
					names.arg=HIST.qdata$mids,cex.names=0.8*cex,
					xlab=pred, ylab="proportion")
			mtext("Training Data", side=3, line=2.5, cex=1.25*cex, font=2 )
			#axis(1,at=seq(min(BP),max(BP),5))
			legend(	x="topright",inset=c(-.3,-.2),title=response.name,bg="white",
					legend=rev(response.colors$category),fill=rev(response.colors$colors),cex=0.6*cex)}
			par(opar)
		}


		### ### ### SCREEN 2 ### ### ###
		#set screen to right half. upper right
		#print("Screen 2 print")		

		opar <- par(cex=cex, mgp=c(1.5, .5, 0),mar=c(3,3,4,1)+.1)
		BXP<-bxp(bp,show.names=FALSE)
		axis(side=3, at=BXP, labels=bp$names,cex.axis=1.1*cex,tick=FALSE,line=0,font=2)
		mtext(pred,side=2,line=1.8,cex=cex)
		par(opar)

		### ### ### SCREEN 3 ### ### ###
		#set screen to right half, lower left
		#print("Screen 3 print")

		opar <- par(cex=cex, mgp=c(1.5, .5, 0),mar=c(3,3,4,1)+.1)
		hist(qdata[,pred],xlim=zlim, xlab=pred, main="Training Data",col=hcol.qdata)
		par(opar)

		### ### ### SCREEN 4 ### ### ###
		#set screen to right half, lower right
		#print("Screen 4 print")

		opar <- par(cex=cex, mgp=c(1.5, .5, 0),mar=c(3,3,4,1)+.1)
		#hist(pred.rast,xlim=zlim, xlab=pred, main="Raster",col=hcol.rast)
		hist(pred.rast,xlim=zlim, xlab=pred, main="Raster",col=hcol.rast)
		par(opar)

		### ### ### SCREEN 5 ### ### ###
		#set screen to left half
		#print("Screen 5 print")

		opar<-par(cex=cex,mar=c(4,3,4,1)+.1)

		image( img.mask, 
		       col = c("grey80","black"),
		       xlab="", ylab="", 
		       zlim=c(0,1),
		       asp=1, bty="n", main="")

		image( img.data, 
		       col = col.ramp,
		       #xlab="", ylab="", 
		       zlim=zlim,
		       #asp=1, bty="n", main="",
			 add=TRUE)


		legend( #x=xmax(mapgrid),y=ymax(mapgrid),
			  x="topright",
			  bg="white",
		        legend=c(legend.label,"NA","MASK"),
		        fill=c(legend.colors,"grey90","black"),
		        #bty="n",
		        cex=cex)
		mtext(paste("Extrapolation Mask -",pred),side=3,line=2.2,cex=1.5*cex)
		mtext("(masked to range of training)",line=1,cex=cex)
		
		mtext(paste("Training Range: ",signif(TRAIN.min,4)," to ",signif(TRAIN.max,4),
				"          Raster Range: ",signif(RAST.min,4)," to ",signif(RAST.max,4),sep=""),
			side=1, line=3,cex=cex,adj=0)

		par(opar)

		if(!device.type[i]%in%c("default","none")){dev.off()}
	}
}
}



#############################################################################################
#############################################################################################
################################# Require Functions #########################################
#############################################################################################
#############################################################################################

###model.type=="CF"###
REQUIRE.party<-function (stopIfAbsent = TRUE) 
{
    y <- getOption("partyLoaded")
    w <- getOption("warn")
    options(warn = -1)
    x <- isTRUE(try(requireNamespace("party", quietly = TRUE)))
    options(warn = w)
    if (!isTRUE(y)) {
        if (x) {
            options(partyLoaded = TRUE)
		#library("party")
            return(TRUE)
        }
        else if (stopIfAbsent) {
            stop("package 'party' used for model type 'CF' is not available")
        }
        else {
            return(FALSE)
        }
    }
    return(TRUE)
}

###model.type=="QRF"###
REQUIRE.quantregForest<-function (stopIfAbsent = TRUE) 
{
    y <- getOption("quantregForestLoaded")
    w <- getOption("warn")
    options(warn = -1)
    x <- isTRUE(try(requireNamespace("quantregForest", quietly = TRUE)))
    options(warn = w)
    if (!isTRUE(y)) {
        if (x) {
            options(quantregForestLoaded = TRUE)
		#library("quantregForest")
            return(TRUE)
        }
        else if (stopIfAbsent) {
            stop("package 'quantregForest' used for model type 'QRF' is not available")
        }
        else {
            return(FALSE)
        }
    }
    return(TRUE)
}

###model.type=="SGB" or "QSGB"###
REQUIRE.gbm<-function (stopIfAbsent = TRUE) 
{
    y <- getOption("gbmLoaded")
    w <- getOption("warn")
    options(warn = -1)
    x <- isTRUE(try(requireNamespace("gbm", quietly = TRUE)))
    options(warn = w)
    if (!isTRUE(y)) {
        if (x) {
            options(gbmLoaded = TRUE)
		#library("gbm")
		#if(!"gbm"%in%loadedNamespaces()){attachNamespace("gbm")}
            return(TRUE)
        }
        else if (stopIfAbsent) {
            stop("package 'gbm' used for model types 'SGB' and 'QSGB' is not available")
        }
        else {
            return(FALSE)
        }
    }
    return(TRUE)
}


#############################################################################################
#############################################################################################
################################## Conditional Forests ######################################
#############################################################################################
#############################################################################################

########## convert vote predictions #################

CF.list2df<-function(pred){
	pred.names<-colnames(pred[[1]])
	pred<-t(sapply(pred,I))
	colnames(pred)<-pred.names
	return(pred)
}

############### convert integer data to numeric #######################

CF.int2num<-function(qdata){
	
	IS.NUM<-sapply(qdata,is.numeric)
	qdata[,IS.NUM]<-sapply(qdata[,IS.NUM],as.numeric)

	return(qdata)
}


#############################################################################################
#############################################################################################
################################## Check functions ##########################################
#############################################################################################
#############################################################################################

########## check device.type #############

check.device.type.nodefault<-function(device.type){

	if(is.null(device.type)){
		device.type <- select.list(c("jpeg","none","pdf","png","postscript","tiff"), title="Diagnostic Output?", multiple = TRUE)
		device.type <- c(device.type,"jpeg")
	}
	if(length(device.type)==0 || is.null(device.type)){
		device.type <- "jpeg"
	}

	if(!is.null(device.type)){
		if(any(device.type%in%c("default","windows"))){
			stop(	"To enable 'default' on screen graphics set 'allow.default.graphics=TRUE'.\n", 
				"Use with caution.\n",
				"If graphics window is closed or moved while R is in process of creating graphic,\n", 
				"it may crash entire R session.")
		}else{
			if(any(!device.type%in%c("jpeg","none","pdf","png","postscript","tiff"))){
				stop("illegal 'device.type' device types must be one or more of 'jpeg' 'pdf' 'png' 'postscript' 'tiff' or 'none'")
			}
			device.type<-sort(device.type)
		}
	}

	if("none"%in%device.type){device.type<-"none"}

	return(device.type)
}


########## check device.type #############

check.device.type<-function(device.type){

	if(is.null(device.type)){
		device.type <- select.list(c("default","jpeg","none","pdf","png","postscript","tiff"), title="Diagnostic Output?", multiple = TRUE)
		device.type <- c(device.type,"default")
	}
	if(length(device.type)==0 || is.null(device.type)){
		device.type <- "default"
	}

	if(!is.null(device.type)){
		device.type[device.type=="windows"]<-"default"
		if(any(!device.type%in%c("default","jpeg","none","pdf","png","postscript","tiff"))){
			stop("illegal 'device.type' device types must be one or more of 'default' 'jpeg' 'pdf' 'png' 'postscript' 'tiff' or 'none'")
		}
		device.type<-sort(device.type)
		if("default"%in%device.type){
			device.type<-c(device.type[device.type!="default"],"default")
		}
	}

	if("none"%in%device.type){device.type<-"none"}

	return(device.type)
}


########## check predList #################

check.predList<-function(model.obj,model.type){
if("QRF"%in%names(model.obj)){
	predList<-row.names(model.obj$QRF$importance)
}else{
	if(model.type == "RF" ){predList<-row.names(model.obj$importance)}
	if(model.type == "QRF"){predList<-row.names(model.obj$importance)}
	if(model.type == "CF"){inputs<-model.obj@data@get("input");predList<-colnames(inputs)}
	if(model.type == "SGB"){predList<-model.obj$var.names}
}
return(predList)
}

########## check model.levels #################

check.model.levels<-function(model.obj,model.type){

model.levels<-NULL

if("QRF"%in%names(model.obj)){
	var.factors<-!sapply(model.obj$QRF$forest$xlevels,identical,0)
		if( any(var.factors)){
			model.levels<-as.list(1:sum(var.factors))
			names(model.levels)<-names(var.factors)[var.factors]
			for(p in 1:sum(var.factors)){
				model.levels[[p]]<-model.obj$QRF$forest$xlevels[var.factors][[p]]}}
}else{
	if(model.type=="SGB"){
		var.factors<-model.obj$var.type!=0
		if( any(var.factors)){
			model.levels<-as.list(1:sum(var.factors))
			names(model.levels)<-model.obj$var.names[var.factors]
			for(p in 1:sum(var.factors)){
				model.levels[[p]]<-model.obj$var.levels[var.factors][[p]]}}}

	if(model.type=="RF"){
		var.factors<-!sapply(model.obj$forest$xlevels,identical,0)
		if( any(var.factors)){
			model.levels<-as.list(1:sum(var.factors))
			names(model.levels)<-names(var.factors)[var.factors]
			for(p in 1:sum(var.factors)){
				model.levels[[p]]<-model.obj$forest$xlevels[var.factors][[p]]}}}

	if(model.type == "CF"){
		inputs<-model.obj@data@get("input")
		var.factors<-sapply(inputs,is.factor)
		if( any(var.factors)){
			model.levels<-as.list(1:sum(var.factors))
			names(model.levels)<-names(var.factors)[var.factors]
			for(p in 1:sum(var.factors)){
				model.levels[[p]]<-levels(inputs[[names(model.levels)[p]]])}}}

	if(model.type=="QRF"){
		var.factors<-!sapply(model.obj$forest$xlevels,identical,0)
		if( any(var.factors)){
			model.levels<-as.list(1:sum(var.factors))
			names(model.levels)<-names(var.factors)[var.factors]
			for(p in 1:sum(var.factors)){
				model.levels[[p]]<-model.obj$forest$xlevels[var.factors][[p]]}}}
}
return(model.levels)
}

########## check model type #################

check.model.type<-function(model.obj){

if("QRF"%in%names(model.obj)){
	model.type<-"QRF"
}else{
	model.type.long<-attr(model.obj,"class")
	if("randomForest"%in%model.type.long){
		if("quantregForest"%in%model.type.long){
			model.type<-"QRF"
		}else{
			model.type<-"RF"
		}
	}else{
		if("RandomForest"%in%model.type.long){
			model.type<-"CF"
		}else{
			if("gbm"%in%model.type.long){
				model.type<-"SGB"
			}else{
				model.type<-"unknown"
			}
		}
	}
}

	if(model.type=="unknown"){
		stop("model.obj is of unknown type")}

	return(model.type)
}

############ check response.type ##################

check.response.type<-function(model.obj,model.type,ONEorTWO){

	if(model.type=="RF"){
		response.type<-switch(model.obj$type,"regression"="continuous","classification"="classification","unknown")
		if(response.type=="classification"){
			if(identical(levels(model.obj$y),c("0","1"))){
				response.type<-"binary"
			}else{
				response.type<-"categorical"}}}
	
	if(model.type=="QRF"){response.type<-"continuous"}
	
	if(model.type=="CF"){    
		if(all(model.obj@responses@is_nominal)){
			if(identical(model.obj@responses@levels[[1]],c("0","1"))){
				response.type<-"binary"
			}else{
				response.type<-"categorical"
			}
		}else{
			response.type<-"continuous"}}
	
	if(model.type=="SGB"){
		response.type<-switch(model.obj$distribution$name,"gaussian"="continuous","bernoulli"="binary","multinomial"="categorical","unknown")}
	
	if(response.type=="unknown"){stop("supplied ", ONEorTWO," has an unknown response type")}

	return(response.type)
}

############ check response.name ##################

check.response.name<-function(model.obj,model.type,response.name,qdata=NULL){

	model.response.name<-NULL
	if(model.type=="CF"){
		if(length(model.obj@responses@variables)==1){
			model.response.name <- names(model.obj@responses@variables)
		}else{
			stop("ModelMap does not support multivariate response")
		}
		#response.name <- as.character(model.obj@data@formula$response[2])
	}else{
		if("QRF"%in%names(model.obj)){
			if(!is.null(model.obj$QRF$response)){model.response.name<-model.obj$QRF$response}
		}else{
			if(!is.null(model.obj$response)){model.response.name<-model.obj$response}
		}
	}

	if(!is.null(model.response.name)){
		if(is.null(response.name)){
			response.name<-model.response.name
		}else{
			if(response.name!=model.response.name){
				stop("supplied model.obj has response ",model.response.name," not supplied response ",response.name)
			}
		}
	}

	## If the response variable is NULL, then the user selects variable from pop-up list.
	if (is.null(response.name)){
		if(!is.null(qdata)){
			if(is.na(qdata)){
				response.name <- "response"
			}else{
				response.name <- select.list(names(qdata), title="Select response name.")
				if(response.name=="" || is.null(response.name)){stop("'response.name' is needed")}
			}	
		}else{stop("'response.name' is needed")}
	}

return(response.name)
}

#############################################################################################
#############################################################################################
################################## Importance Plot ##########################################
#############################################################################################
#############################################################################################

########## extract importance SGB ############

imp.extract.sgb<-function(model.obj,imp.type=NULL){
	imp.type.gbm<-switch(imp.type,"1"=gbm::permutation.test.gbm,"2"=gbm::relative.influence)
	IMP.SGB<-gbm::summary.gbm(model.obj, method=imp.type.gbm, plotit=FALSE)
	names(IMP.SGB)<-c("pred","imp")
	row.names(IMP.SGB)<-IMP.SGB$pred
	IMP.SGB<-IMP.SGB[order(IMP.SGB$imp,decreasing=FALSE),]
	return(IMP.SGB)
}

######### extract importance RF ############

imp.extract.rf<-function(model.obj,imp.type=NULL,class=NULL){
	IMP.RF<-importance(model.obj,type=imp.type,class=class)
	IMP.RF<-data.frame(pred=rownames(IMP.RF),imp=IMP.RF[,1])
	IMP.RF<-IMP.RF[order(IMP.RF[,2],decreasing=FALSE),]
	return(IMP.RF)
}

######### extract importance QRF ############

imp.extract.qrf<-function(model.obj,imp.type=1,ONEorTWO="model.obj"){

if("QRF"%in%names(model.obj)){
	if(!"quantiles"%in%names(model.obj$QRF)){stop("QRF model ",ONEorTWO," was built with 'importance=FALSE'")}
	IMP.RF<-importance(model.obj$RF)
	IMP.RF<-data.frame(pred=rownames(IMP.RF),RF=IMP.RF)
	IMP.RF<-IMP.RF[order(IMP.RF[,imp.type],decreasing=FALSE),]

	IMP.QRF<-importance(model.obj$QRF)
	IMP.QRF<-IMP.QRF[match(IMP.RF$pred,row.names(IMP.QRF)),]

	IMP.QRF<-cbind(IMP.RF,IMP.QRF)

}else{
	if(!"quantiles"%in%names(model.obj)){stop("QRF model ",ONEorTWO," was built with 'importance=FALSE'")}
	IMP.QRF<-importance(model.obj)
	IMP.QRF<-cbind(data.frame(pred=rownames(IMP.QRF)),IMP.QRF)
	IMP.QRF<-IMP.QRF[order(IMP.QRF[,2],decreasing=FALSE),]

}
return(IMP.QRF)
}

######### extract importance CF ############

imp.extract.cf<-function(	model.obj,
					imp.type=NULL,
					mincriterion = 0, 
					conditional = FALSE, 
       				threshold = 0.2, 
					nperm = 1){
	if(imp.type==1){
		IMP.CF<-party::varimp(model.obj,mincriterion=mincriterion,conditional=conditional,threshold=threshold,nperm=nperm)}
	if(imp.type==2){
		IMP.CF<-party::varimpAUC(model.obj,mincriterion=mincriterion,conditional=conditional,threshold=threshold,nperm=nperm)}
	IMP.CF<-data.frame(pred=names(IMP.CF),imp=IMP.CF)
	IMP.CF<-IMP.CF[order(IMP.CF[,2],decreasing=FALSE),]
	return(IMP.CF)
}

######### scale importance ###############

imp.scale<-function(IMP,scale.by="max"){

	if(scale.by=="max"){
		IMP$imp[IMP$imp<0]<-0
		IMP$imp<-(IMP$imp/max(IMP$imp))*1
	}

	if(scale.by=="sum"){
		IMP$imp[IMP$imp<0]<-0
		IMP$imp<-(IMP$imp/sum(IMP$imp))*1
	}

	return(IMP)
}



#############################################################################################
#############################################################################################
###################################### Get Rasts ############################################
#############################################################################################
#############################################################################################

getRasts<-function(){
	## This function prompts user to browse to and select raster layers used for model predictors.
	## Returns vector of rasters.  

	## Adds to file filters to Cran R Filters table.
	Filters<-rbind(Filters,img=c("Imagine files (*.img)", "*.img"))
	Filters<-rbind(Filters,csv=c("Comma-delimited files (*.csv)", "*.csv"))

	userPrompt <- "Yes"
      rastnmVector <- {}
      while(userPrompt == "Yes"){
		## Gets raster format from user.
      	rasterType = select.list(c("Imagine Image", "ArcInfo Grid"), title="Select type of raster.")
		if(rasterType==""){stop("Type of Raster must be selected")}
		if(rasterType == "Imagine Image"){	
			rasts <- choose.files(caption="Select image", filters = Filters["img",], multi = TRUE)
				if(is.null(rasts)){
					stop("")}
			}
		if(rasterType == "ArcInfo Grid"){
			rasts = choose.dir(default=getwd(),caption="Select grid")
			if(is.null(rasts)){
				stop("")}
		}
		if(is.null(rasterType)){
			stop("")
		}

		## Prompts user to select more.
		userPrompt = select.list(c("Yes", "No"), title="Another?")
		if(userPrompt==""){userPrompt<-"No"}
		## Compiles list of rasters.
		rastnmVector = c(rastnmVector, rasts)
	}
	return(rastnmVector)
}

#############################################################################################
#############################################################################################
###################################### Diagnostics ##########################################
#############################################################################################
#############################################################################################

diagnostics.function<-function(	model.obj=NULL,
						model.type,
						predList,
						predFactor,
						PRED,
						qdata,
						MODELfn=NULL,
						MODELpredfn=NULL,
						main=basename(MODELpredfn),
						response.name,
						response.type,
						prediction.type,
						quantiles=NULL,
						imp.quantiles=NULL,
						LOWER=0.10,
						UPPER=0.90,
						folder=getwd(),
						device.type="default",
						res=72,
						device.width=7,
						device.height=7,
						units="in",
						pointsize=pointsize,
						cex=par$cex,
						req.sens,
						req.spec,
						FPC,
						FNC	){

#warning("diagnostic warning")
### Write out tables ###

if(response.type == "binary"){
	OPTTHRESHfn<-paste(MODELpredfn,"_optthresholds.csv",sep="")
	PREDPREVfn<-paste(MODELpredfn,"_prevalence.csv",sep="")

	opt.thresh<-suppressWarnings(error.threshold.plot(	PRED,opt.methods=optimal.thresholds(),plot.it=FALSE,
										req.sens=req.sens,req.spec=req.spec,FPC=FPC,FNC=FNC))

	pred.prev<-predicted.prevalence(PRED, threshold = opt.thresh$threshold)
	pred.prev<-cbind(opt.thresh$opt.methods, pred.prev)

	write.table(	opt.thresh,file=OPTTHRESHfn,sep=",",row.names=FALSE)
	write.table(	pred.prev,file=PREDPREVfn,sep=",",row.names=FALSE)
}

if(response.type == "categorical"){

	CMXfn<-paste(MODELpredfn,"_cmx.csv",sep="")

	LEVELS.pred <- levels(PRED$pred)
	LEVELS.obs  <- levels(PRED$obs)

	print(paste("LEVELS.pred: ",paste(LEVELS.pred, collapse=" ")))
	print(paste("LEVELS.obs:  ",paste(LEVELS.obs,  collapse=" ")))


	LEVELS<-sort(unique(c(LEVELS.pred,LEVELS.obs)))

	CMX<-table(	predicted = factor(PRED$pred,levels=LEVELS),
			observed = factor(PRED$obs,levels=LEVELS))


	CMX.out<-matrix("cmx",nrow(CMX)+5,ncol(CMX)+4)
	CMX.out[1,(1:ncol(CMX))+2]<-"observed"
	CMX.out[2,(1:ncol(CMX))+2]<-colnames(CMX)

	CMX.out[(1:nrow(CMX))+2,1]<-"predicted"
	CMX.out[(1:nrow(CMX))+2,2]<-rownames(CMX)

	CMX.out[(1:nrow(CMX))+2,(1:ncol(CMX))+2]<-CMX

	###Kappa###
	KAPPA<-signif(Kappa(CMX),6)
	CMX.out[1,1:2]<-names(KAPPA)
	CMX.out[2,1]<-KAPPA[1,1]
	CMX.out[2,2]<-KAPPA[1,2]

	###Totals###
	CMX.out[1:2,ncol(CMX.out)-1]<-"total"
	CMX.out[nrow(CMX.out)-2,1:2]<-"total"
	
	CMX.out[(1:nrow(CMX))+2,ncol(CMX.out)-1]<-apply(CMX,1,sum)
	CMX.out[nrow(CMX.out)-2,(1:ncol(CMX))+2]<-apply(CMX,2,sum)
	CMX.out[nrow(CMX.out)-2,ncol(CMX.out)-1]<-sum(CMX)

	###marginals###
	CMX.diag<-diag(CMX)

	CMX.out[1:2,ncol(CMX.out)]<-"Commission"
	CMX.out[nrow(CMX.out)-1,1:2]<-"Omission"

	CMX.out[(1:nrow(CMX))+2,ncol(CMX.out)]<-1-(CMX.diag/apply(CMX,1,sum))
	CMX.out[nrow(CMX.out)-1,(1:ncol(CMX))+2]<-1-(CMX.diag/apply(CMX,2,sum))

	###pcc###
	CMX.out[nrow(CMX.out)-2,ncol(CMX.out)]<-"PCC"
	CMX.out[nrow(CMX.out)-1,ncol(CMX.out)-1]<-"PCC"

	CMX.out[nrow(CMX.out)-1,ncol(CMX.out)]<-sum(CMX.diag)/sum(CMX)

	###MAUC###

	#if(prediction.type=="CV"){
      #      PRED.mauc = PRED[4:(ncol(PRED)-1)]
	#}else{
	#	PRED.mauc = PRED[,4:ncol(PRED)]
	#}

	#PRED.mauc <- PRED[,LEVELS.pred]

	if(prediction.type=="CV"){
		obs.in.pred <- LEVELS.obs[LEVELS.obs%in%names(PRED)[4:(ncol(PRED)-1)]]
	}else{
		obs.in.pred <- LEVELS.obs[LEVELS.obs%in%names(PRED)[4:ncol(PRED)]]}
	obs.in.pred <- obs.in.pred[apply(PRED[,obs.in.pred],2,sum)!=0]

	PRED.mauc <- PRED[,obs.in.pred]

	#if(any(!LEVELS.obs%in%LEVELS.pred)){
	if(any(!LEVELS.obs%in%obs.in.pred)){

		LEVELS.new <- LEVELS.obs[!LEVELS.obs%in%LEVELS.pred]
		LEVELS.new.paste <-paste(LEVELS.new,collapse=" ")

		warning("Response categories: ", LEVELS.new.paste, " in observed data have 0% probability for all datapoints")
		
		PRED.new  <- matrix(NA,nrow=nrow(PRED.mauc),ncol=length(LEVELS.new))
		colnames(PRED.new)<-LEVELS.new

		PRED.mauc<-cbind(PRED.mauc,PRED.new)
	}

	VOTE <- multcap(  response = PRED$obs,
                		predicted= as.matrix(PRED.mauc) )

	MAUC  <- HandTill2001::auc(VOTE)

	CMX.out[nrow(CMX.out),1]<-"MAUC"
	CMX.out[nrow(CMX.out),2]<-MAUC

	write.table(CMX.out,file=CMXfn,sep=",",row.names=FALSE,col.names = FALSE)

}


if(response.type == "continuous" && model.type!="QRF"){

	COR.pearson<-cor(PRED$pred,PRED$obs,method="pearson")
	COR.pearson<-round(COR.pearson,2)
	
	COR.spearman<-cor(PRED$pred,PRED$obs,method="spearman")
	COR.spearman<-round(COR.spearman,2)

	Resid<-PRED$obs-PRED$pred
	n<-length(Resid)
	MSE<-mean(Resid^2)
	RMSD<-(sum((Resid^2))/(n-1))^.5
	RMSD<-round(RMSD,2)

	lm.pred<-lm(obs~pred,data=PRED)
	b<-round(lm.pred$coefficients[1],2)
	m<-round(lm.pred$coefficients[2],2)
}

if(model.type=="QRF" && "RFmean"%in%names(PRED)){

	COR.pearson<-cor(PRED$RFmean,PRED$obs,method="pearson")
	COR.pearson<-round(COR.pearson,2)
	
	COR.spearman<-cor(PRED$RFmean,PRED$obs,method="spearman")
	COR.spearman<-round(COR.spearman,2)

	Resid<-PRED$obs-PRED$RFmean
	n<-length(Resid)
	MSE<-mean(Resid^2)
	RMSD<-(sum((Resid^2))/(n-1))^.5
	RMSD<-round(RMSD,2)

	lm.pred<-lm(obs~RFmean,data=PRED)
	b<-round(lm.pred$coefficients[1],2)
	m<-round(lm.pred$coefficients[2],2)
}

### make graphs ###

#
if(!"none"%in%device.type){
for(i in 1:length(device.type)){
#print(paste("Device Type:",device.type[i]))

	main=basename(MODELpredfn)

	### Output filenames ###

	SCATTERPLOTfn<-paste(MODELpredfn,"_scatterplot",sep="")
	IMPORTANCEfn<-paste(MODELpredfn,"_importance",sep="")
	ERRORfn<-paste(MODELpredfn,"_error",sep="")
	THRESHOLDPLOTSfn<-paste(MODELpredfn,"_thresholdplots",sep="")

	SUMMARYfn<-paste(MODELpredfn,"_quantreg",sep="")

	QUANTWHISKfn<-paste(MODELpredfn,"_quantreg_whisker",sep="")
	QUANTCONfn<-paste(MODELpredfn,"_quantreg_confidence",sep="")
	QUANTCONLOGfn<-paste(MODELpredfn,"_quantreg_confidence_log",sep="")

	###################   model.obj Based   ##############################


	if(model.type=="RF"){

		
		###Importance Plot###
		#print(paste("IMPORTANCEfn =",IMPORTANCEfn)

		initialize.device(	PLOTfn=IMPORTANCEfn,DEVICE.TYPE=device.type[i],
						res=res,device.width=device.width,device.height=device.height,
						units=units,pointsize=pointsize,cex=cex)

		opar<-par(cex=cex)
		varImpPlot(model.obj,main="Relative Influence",cex=cex)
		mtext(main,side=3,line=-4,cex=1.3*cex,outer=TRUE)
		par(opar)

		if(!device.type[i]%in%c("default","none")){dev.off()}

		###Importance Plot - category graphs###
		#print("     Importance plots")
		if(response.type%in%c("binary","categorical")){
			Nplots<-table(model.obj$y)

			CATnames<-names(Nplots)
			CATfn<-paste("Category",CATnames,sep="_")
			CATfn<-make.names(CATfn)
			Nfig<-length(CATnames)

			for(j in 1:Nfig){

				IMPFIGfn<-paste(IMPORTANCEfn,CATfn[j],sep="_")

				initialize.device(	PLOTfn=IMPFIGfn,DEVICE.TYPE=device.type[i],
								res=res,device.width=device.width,device.height=device.height,
								units=units,pointsize=pointsize,cex=cex)

				opar<-par(cex=cex)			

				varImpPlot(model.obj,cex=cex,type=1,class=CATnames[j],main="")
				mtext(main,side=3,line=-2,cex=1.3*cex,outer=TRUE)
				mtext(paste("Relative Influence -",CATnames[j],"-",Nplots[j],"plots"),side=3,line=-3.5,cex=1.3*cex,outer=TRUE)

				par(opar)

				if(!device.type[i]%in%c("default","none")){dev.off()}
			}
		}


		###Error Plot###
		#print(     "Error Plots")
		if(response.type == "continuous"){
			initialize.device(	PLOTfn=ERRORfn,DEVICE.TYPE=device.type[i],
							res=res,device.width=device.width,device.height=device.height,
							units=units,pointsize=pointsize,cex=cex)
			opar<-par(cex=cex)
		
			Nmax<-length(model.obj$mse)
			Nplots<-length(model.obj$predicted)

			plot(	1:Nmax,
				model.obj$mse,
				type="l", 
				xlab="ntree",ylab="MSE")
			mtext(main,side=3,line=-2,cex=1.3*cex,outer=TRUE)
			mtext(paste("OOB -",Nplots,"plots"),side=3,line=-3.5,cex=1.1*cex,outer=TRUE)

			par(opar)
			if(!device.type[i]%in%c("default","none")){dev.off()}
		}

		###Error Plot - category graphs###
		if(response.type%in%c("binary","categorical")){
			Nplots<-table(model.obj$y)
			Nplots<-c(sum(Nplots),Nplots)

			CATnames<-colnames(model.obj$err.rate)
			CATfn<-CATnames
			CATfn[-1]<-paste("Category",CATnames[-1],sep="_")
			CATfn<-make.names(CATfn)
			Nfig<-length(CATnames)
			Nmax<-nrow(model.obj$err.rate)

			for(j in 1:Nfig){

				ERRORFIGfn<-paste(ERRORfn,CATfn[j],sep="_")
				initialize.device(	PLOTfn=ERRORfn,DEVICE.TYPE=device.type[i],
								res=res,device.width=device.width,device.height=device.height,
								units=units,pointsize=pointsize,cex=cex)
								opar<-par(cex=cex)
				plot(	1:Nmax,
					model.obj$err.rate[,j],
					type="l", 
					xlab="ntree",ylab="OOB err.rate")
				mtext(main,side=3,line=-2,cex=1.3*cex,outer=TRUE)
				mtext(paste("OOB -",CATnames[j],"-",Nplots[j],"plots"),side=3,line=-3.5,cex=1.1*cex,outer=TRUE)

				par(opar)
				if(!device.type[i]%in%c("default","none")){dev.off()}
			}
		}
	}

	if(model.type=="QRF"){
		###Summary Plot###
		#print(     "Summary Plot")
		initialize.device(	PLOTfn=SUMMARYfn,DEVICE.TYPE=device.type[i],
						res=res,device.width=device.width,device.height=device.height,
						units=units,pointsize=pointsize,cex=cex)
		opar<-par(cex=cex)
		if("QRF"%in%names(model.obj)){plot(model.obj$QRF)}else{plot(model.obj)}
		par(opar)
		if(!device.type[i]%in%c("default","none")){dev.off()}

		###QRF Importance plot###
		if(!is.null(imp.quantiles)){
 			initialize.device(	PLOTfn=paste(IMPORTANCEfn,"_QRF",sep=""),DEVICE.TYPE=device.type[i],
							res=res,device.width=device.width,device.height=device.height,
							units=units,pointsize=pointsize,cex=cex)

			opar<-par(cex=cex)
			if("QRF"%in%names(model.obj)){
				quantregForest::varImpPlot.qrf(model.obj$QRF,main="Quantile Importance",cex=cex)
			}else{
				quantregForest::varImpPlot.qrf(model.obj,main=paste(main,"Quantile Importance"),cex=cex)
			}
			
			mtext(main,side=3,line=.2,cex=1.3*cex)
			par(opar)
			if(!device.type[i]%in%c("default","none")){dev.off()}
		}
		###RF Importance plot###
		if("RF"%in%names(model.obj)){

			initialize.device(	PLOTfn=paste(IMPORTANCEfn,"_RF",sep=""),DEVICE.TYPE=device.type[i],
							res=res,device.width=device.width,device.height=device.height,
							units=units,pointsize=pointsize,cex=cex)

			opar<-par(cex=cex)
			varImpPlot(model.obj$RF,main="Random Forest Model",cex=cex)

			mtext(main,side=3,line=-4,cex=1.3*cex,outer=TRUE)
			par(opar)

			if(!device.type[i]%in%c("default","none")){dev.off()}
		}
	}

	if(model.type=="SGB"){
	
		initialize.device(	PLOTfn=IMPORTANCEfn,DEVICE.TYPE=device.type[i],
						res=res,device.width=device.width,device.height=device.height,
						units=units,pointsize=pointsize,cex=cex)

		opar<-par(las=1,mar=(c(5, 11, 4, 2) + 0.1),cex=cex)
		gbm::summary.gbm(model.obj)
		par(las=0)
		mtext("Relative Influence",side=3,line=.7,cex=1.5*cex)
		mtext(main,side=3,line=2.7,cex=1.5*cex)
		mtext("Predictors",side=2,line=10,cex=1*cex)
		par(opar)

		if(!device.type[i]%in%c("default","none")){dev.off()}
	}


	###################   PRED Based   ##############################

	### binary ###

	if(response.type == "binary"){

		initialize.device(	PLOTfn=THRESHOLDPLOTSfn,DEVICE.TYPE=device.type[i],
						res=res,device.width=device.width,device.height=device.height,
						units=units,pointsize=pointsize,cex=cex)

		opar<-par(cex=cex)
		suppressWarnings(presence.absence.summary(PRED,main=main,legend.cex=cex,opt.legend.cex=cex))
		par(opar)

		if(!device.type[i]%in%c("default","none")){dev.off()}	
	}

### Continuous ###

	if(response.type == "continuous" && model.type!="QRF"){

		initialize.device(	PLOTfn=SCATTERPLOTfn,DEVICE.TYPE=device.type[i],
						res=res,device.width=device.width,device.height=device.height,
						units=units,pointsize=pointsize,cex=cex)


		opar<-par(pty="s",cex=cex)
		lim<-range(PRED$obs,PRED$pred)
		plot(PRED$pred,PRED$obs,xlab="predicted",ylab="observed",xlim=lim,ylim=lim,main="")
		abline(a=0,b=1,lty=2)
		abline(lm.pred)

		mtext(main,side=3,line=1.5,cex=1.3*cex)

		mtext(paste("RMSD:",RMSD," ",sep=""),side=1,line=-4.5,adj=1,cex=.8*cex)	
		mtext(paste("pearson's cor: ",COR.pearson," ",sep=""),side=1,line=-3.5,adj=1,cex=.8*cex)
		mtext(paste("spearman's cor: ",COR.spearman," ",sep=""),side=1,line=-2.5,adj=1,cex=.8*cex)
		mtext(paste("obs = ",m,"(pred) + ",b," ",sep=""),side=1,line=-1.5,adj=1,cex=.8*cex)

		par(opar)

		if(!device.type[i]%in%c("default","none")){dev.off()}
	}	
	if(model.type=="QRF"){
		main.quant<-paste( round((UPPER-LOWER)*100,0),"% Confidence Interval",sep="")

		###Scatter Plot###

		initialize.device(	PLOTfn=SCATTERPLOTfn,DEVICE.TYPE=device.type[i],
						res=res,device.width=device.width,device.height=device.height,
						units=units,pointsize=pointsize,cex=cex)

		opar<-par(pty="s",cex=cex)
		if("RFmean"%in%names(PRED)){lim<-range(PRED$obs,PRED$RFmean,PRED$median)}else{range(PRED$obs,PRED$median)}
		plot(PRED$median,PRED$obs,xlab="predicted",ylab="observed",xlim=lim,ylim=lim,main="",col="blue",pch=1)
		abline(a=0,b=1,lty=2)
			
		mtext(main,side=3,line=1.5,cex=1.3*cex)
			
		if("RFmean"%in%names(PRED)){
			points(PRED$RFmean,PRED$obs,col="red",pch=1)
			abline(lm.pred)
			mtext("RF model:",side=1,line=-5.5,adj=1,cex=.8*cex,font=2)
			mtext(paste("RMSD:",RMSD," ",sep=""),side=1,line=-4.5,adj=1,cex=.8*cex)	
			mtext(paste("pearson's cor: ",COR.pearson," ",sep=""),side=1,line=-3.5,adj=1,cex=.8*cex)
			mtext(paste("spearman's cor: ",COR.spearman," ",sep=""),side=1,line=-2.5,adj=1,cex=.8*cex)
			mtext(paste("obs = ",m,"(pred) + ",b," ",sep=""),side=1,line=-1.5,adj=1,cex=.8*cex)
			legend("topright",legend=c("QRF median","RF mean"),col=c("blue","red"),pch=1,cex=.8*cex,inset=0.02,bg="white")
		}else{
			legend("topright",legend=c("QRF median"),col=c("blue"),pch=1,cex=.8*cex,inset=0.02,bg="white")
		}
		par(opar)
	
		if(!device.type[i]%in%c("default","none")){dev.off()}

		###Whisker Plot###
		PRED.s<-PRED[order(PRED$obs),]

		initialize.device(	PLOTfn=QUANTWHISKfn,DEVICE.TYPE=device.type[i],
						res=res,device.width=device.width*2,device.height=device.height,
						units=units,pointsize=pointsize,cex=cex)

		#whisker.quantiles<-c(range(quantiles))
		#quantlab<-sprintf("P%02d", 100*whisker.quantiles)
		#main<-paste(as.character(100*(whisker.quantiles[2]-whisker.quantiles[1])),"% Confidence Interval",sep="")
		#print(quantlab)

		opar<- par(mar=par()$mar+c(0,0,0,5))
		lim<-range(PRED.s[,-1])
		plot(	PRED.s$obs,PRED.s$obs,pch=16,
			xlim=lim,ylim=lim,xlab="observed response",ylab="response",main=main.quant)
		for(j in 1:nrow(PRED.s)){
			lines(x=c(PRED.s$obs[j],PRED.s$obs[j]),
				y=c(PRED.s$lower[j],PRED.s$upper[j]),
				col="blue")
		}
		points(PRED.s$obs,PRED.s$obs,pch=16,col="black")
		points(PRED.s$obs,PRED.s$median,pch=16,col="blue")
		if("RFmean"%in%names(PRED)){points(PRED.s$obs,PRED.s$RFmean,pch=16,col="red")}
		op.leg<-par(xpd=NA)
		if("RFmean"%in%names(PRED)){
			legend(	x=lim[2],y=lim[2],
					c("observed","QRF median","RF mean"),
					col=c("black","blue","red"),
					pch=16,bg="white")
		}else{
			legend(	x=lim[2],y=lim[2],
					c("observed","QRF median"),
					col=c("black","blue"),
					pch=16,bg="white")
		}
		par(op.leg)
		par(opar)
		if(!device.type[i]%in%c("default","none")){dev.off()}

		###Predictor Whisker Plots###
		predCont<-predList[!predList%in%predFactor]
		#print(paste("predCont:",predCont))

		for(xvar in predCont){
			#print(paste("      xvar:",xvar))

			ORD<-order(qdata[,xvar])
			PRED.s<-PRED[ORD,]
			qdata.s<-qdata[ORD,]

			initialize.device(	PLOTfn=paste(QUANTWHISKfn,"_",xvar,sep=""),DEVICE.TYPE=device.type[i],
							res=res,device.width=device.width*2,device.height=device.height,
							units=units,pointsize=pointsize,cex=cex)

			#whisker.quantiles<-c(range(quantiles))
			#quantlab<-sprintf("P%02d", 100*whisker.quantiles)
			#main<-paste(as.character(100*(whisker.quantiles[2]-whisker.quantiles[1])),"% Confidence Interval",sep="")
			#print(quantlab)

			opar<- par(mar=par()$mar+c(0,0,0,5))
			xlim<-range(qdata.s[,xvar])
			ylim<-range(PRED.s[,-1],na.rm=TRUE)

			main.whisker=paste(response.name,"~",xvar," (",main.quant,")",sep="")

			plot(	qdata.s[,xvar],PRED.s$obs,pch=16,
				xlim=xlim,ylim=ylim,xlab=xvar,ylab="response",type="n")
			mtext(main.whisker,side=3, line=2, cex=1.5, font=2)
			for(j in 1:nrow(qdata.s)){
				lines(x=c(qdata.s[,xvar][j],qdata.s[,xvar][j]),
					y=c(PRED.s[,"lower"][j],PRED.s[,"upper"][j]),
					col="blue")
			}

			points(qdata.s[,xvar],PRED.s$obs,pch=16,col="black")
			points(qdata.s[,xvar],PRED.s$median,pch=16,col="blue")
			if("RFmean"%in%names(PRED)){points(qdata.s[,xvar],PRED.s$RFmean,pch=16,col="red")}
			op.leg<-par(xpd=NA)
			if("RFmean"%in%names(PRED)){
				legend(	x=xlim[2],y=ylim[2],
						c("observed","QRF median","RF mean"),
						col=c("black","blue","red"),
						pch=16,bg="white")
			}else{
				legend(	x=xlim[2],y=ylim[2],
						c("observed","QRF median"),
						col=c("black","blue"),
						pch=16,bg="white")
			}
			par(op.leg)
			par(opar)
			if(!device.type[i]%in%c("default","none")){dev.off()}
		}

#		###Confidence Interval Plot###
#		PRED.s<-PRED[order(PRED$obs),]
#		PRED.log<-PRED.s
#		PRED.log[,-1][PRED.log[,-1]<=0]<-NA
#		PRED.log[,-1]<-log(PRED.log[,-1],base=10)
#		PRED.log[,-1][is.na(PRED.log[,-1])|PRED.log[,-1]<0]<-0
#	
#		initialize.device(	PLOTfn=QUANTCONfn,DEVICE.TYPE=device.type[i],
#						res=res,device.width=device.width,device.height=device.height,
#						units=units,pointsize=pointsize,cex=cex)
#
#		opar<- par(mar=par()$mar+c(0,0,0,5),pty="s")
#		lim<-range(PRED.s[,-1])
#		plot(	PRED.s$obs,PRED.s$obs,pch=16,
#			xlim=lim,ylim=lim,xlab="observed response",ylab="response",
#			main="90% Confidence Interval")
#		polygon(	c(PRED.s$obs, rev(PRED.s$obs)), 
#				c(PRED.s$lower, rev(PRED.s$upper)), 
#				col = "lightgrey", border = "blue")
#		points(PRED.s$obs,PRED.s$obs,pch=16,col="black")
#		points(PRED.s$obs,PRED.s$median,pch=16,col="blue")
#		if("RFmean"%in%names(PRED)){points(PRED.s$obs,PRED.s$RFmean,pch=16,col="red")}
#		op.leg<-par(xpd=NA)
#		if("RFmean"%in%names(PRED)){
#			legend(	x=lim[2],y=lim[2],
#					c("observed","QRF median","RF mean"),
#					col=c("black","blue","red"),
#					pch=16,bg="white")
#		}else{
#			legend(	x=lim[2],y=lim[2],
#					c("observed","QRF median"),
#					col=c("black","blue"),
#					pch=16,bg="white")
#		}
#		par(op.leg)
#		par(opar)
#		if(!device.type[i]%in%c("default","none")){dev.off()}
#
#
#		###Log transformed Confidence Interval###
#	
#		initialize.device(	PLOTfn=QUANTCONLOGfn,DEVICE.TYPE=device.type[i],
#						res=res,device.width=device.width,device.height=device.height,
#						units=units,pointsize=pointsize,cex=cex)
#
#		opar<- par(mar=par()$mar+c(0,0,0,5),pty="s")
#		xlim<-range(PRED.s[,"obs"])
#		ylim<-range(PRED.log[,-1])
#		plot(	PRED.s$obs,PRED.log$obs,pch=16,
#			xlim=xlim,ylim=ylim,xlab="observed response",ylab="log10(response)",
#			main="Log Transformed 90% Confidence Interval")
#		polygon(	c(PRED.s$obs, rev(PRED.s$obs)), 
#				c(PRED.log$lower, rev(PRED.log$upper)), 
#				col = "lightgrey", border = "blue")
#		points(PRED.s$obs,PRED.log$obs,pch=16,col="black")
#		points(PRED.s$obs,PRED.log$median,pch=16,col="blue")
#		if("RFmean"%in%names(PRED)){points(PRED.s$obs,PRED.log$RFmean,pch=16,col="red")}
#		op.leg<-par(xpd=NA)
#		if("RFmean"%in%names(PRED)){
#			legend(	x=xlim[2],y=ylim[2],
#					col=c("black","blue","red"),
#					pch=16,bg="white")
#		}else{
#			legend(	x=xlim[2],y=ylim[2],
#					c("observed","QRF median"),
#					col=c("black","blue"),
#					pch=16,bg="white")
#		}
#
#		par(op.leg)
#		par(opar)
#		if(!device.type[i]%in%c("default","none")){dev.off()}
##########################################################
	}
}

}
}



#############################################################################################
################################ SGB - Model Creation #######################################
#############################################################################################

model.SGB<-function(	qdata,
				predList,
				response.name,
				response.type,
				seed=NULL,
				n.trees=NULL,                 # number of trees
				shrinkage=0.001,   	      # shrinkage or learning rate,
                 	 	interaction.depth=10,		# 1: additive model, 2: two-way interactions, etc.
				bag.fraction = 0.5,          	# subsampling fraction, 0.5 is probably best
				nTrain = NULL,
				#train.fraction = NULL,       	# fraction of data for training,
                 	 	n.minobsinnode = 10,         	# minimum total weight needed in each node
				keep.data=TRUE,
				var.monotone = NULL
){


## This function generates a model using gbm.
##	Inputs: Full dataset, training indices, predictor names, response name, and seed (optional)
##	Output: SGB model

if(!is.null(seed)){
	set.seed(seed)}

if(response.type=="binary"){distribution="bernoulli"}
if(response.type=="continuous"){distribution="gaussian"}
if(response.type=="categorical"){distribution="multinomial"}

flag.nt<-FALSE
if(response.type=="categorical"){
	flag.nt<-TRUE
	if(is.null(n.trees)){n.trees<-5000}}

qdata.x<-qdata[,match(predList,names(qdata))]

qdata.y<-qdata[,response.name]
if(response.type=="binary"){qdata.y[qdata.y>0]<-1}

if(is.null(n.trees)){

	if(keep.data==FALSE){
		warning("keep.data reset to TRUE because data needed for 'gbm.more()' function needed for OOB determination of optimal number of trees")}

	SGB <- gbm::gbm.fit(	x=qdata.x,
				y=qdata.y,        
				distribution=distribution,
				n.trees=100,                	
				shrinkage=shrinkage, 
				interaction.depth=interaction.depth,		
				bag.fraction = bag.fraction,   
				nTrain = nTrain,       	
				#train.fraction = train.fraction,       	           		
				n.minobsinnode = n.minobsinnode,
				keep.data=TRUE,
				var.monotone=var.monotone)

	# check performance using an out-of-bag estimator
	best.iter <- suppressWarnings(gbm::gbm.perf(SGB,method="OOB",plot.it=FALSE))
      
	# iterate until a sufficient number of trees are fit

	while(SGB$n.trees - best.iter < 10){
     	 	# do 100 more iterations
      	SGB <- gbm::gbm.more(SGB,100)          
      	best.iter <- suppressWarnings(gbm::gbm.perf(SGB,method="OOB",plot.it=FALSE))
	}
	SGB$best.iter <- best.iter
	warning("ModelMap currently uses OOB estimation to determine optimal number of trees in SGB model when calling 'gbm.perf' in the gbm package however OOB generally underestimates the optimal number of iterations although predictive performance is reasonably competitive however using cv.folds>0 when calling gbm usually results in improved predictive performance but is not yet supported in ModelMap")

}else{

	SGB <- gbm::gbm.fit(	x=qdata.x,
				y=qdata.y,        
				distribution=distribution,
				n.trees=n.trees,                	
				shrinkage=shrinkage, 
				interaction.depth=interaction.depth,		
				bag.fraction = bag.fraction, 
				nTrain = nTrain,         	
				#train.fraction = train.fraction,       	           		
				n.minobsinnode = n.minobsinnode,
				keep.data=keep.data,
				var.monotone=var.monotone)

	if(!is.null(nTrain) && nTrain<nrow(qdata.x)){
		SGB$best.iter <- suppressWarnings(gbm::gbm.perf(SGB,method="test",plot.it=FALSE))
		if(SGB$best.iter>0.9*n.trees){
			warning("best number of trees is ", SGB$best.iter, " and total number trees tested was ", n.trees, " therefore you may want to explore increasing the 'n.trees' argument")
		}
	}else{
		if(flag.nt){
			SGB$best.iter <- suppressWarnings(gbm::gbm.perf(SGB,method="OOB",plot.it=FALSE))
			if(SGB$best.iter>0.9*n.trees){
				warning("best number of trees is ", SGB$best.iter, " and total number trees tested was ", n.trees, " therefore you may want to explore increasing the 'n.trees' argument")
			}
		}
	}
}

SGB$response<-response.name

return(SGB)

}




#############################################################################################
################################### SGB - Predict ###########################################
#############################################################################################

prediction.SGB<-function(	prediction.type,
					qdata,
					response.name=deparse(substitute(SGB$response.name)),
					SGB,
					n.trees
					){

## This function makes predictions for SGB model.
##	Inputs: Training data, training data indices, predictor names, response variable name,
##			the Random Forest model and what to do if NAs are in predictors (default).
##	Output: Observed and predicted values.

response.type<-switch(SGB$distribution$name,"gaussian"="continuous","bernoulli"="binary","multinomial"="categorical","unknown")

if(response.type=="unknown"){
	stop("supplied model.obj has an unknown response type")
}

if(is.null(response.name)){
	stop("must provide response name")}

predList<-SGB$var.names

qdata.x<-qdata[,match(predList,names(qdata))]
qdata.y<-qdata[,response.name]

if(response.type=="binary"){qdata.y[qdata.y>0]<-1}

if(prediction.type=="TRAIN"){

	if(!is.null(SGB$levels)){
		for(p in names(SGB$levels)){
			qdata.x[,p]<-factor(qdata.x[,p],levels=SGB$levels[[p]])
		}
	}

	pred<-gbm::predict.gbm(	object=SGB,
					newdata=qdata.x,
					n.trees=n.trees,
					type="response",
					single.tree=FALSE)
}

if(prediction.type=="TEST"){
	
	pred<-gbm::predict.gbm(	object=SGB,
					newdata=qdata.x,
					n.trees=n.trees,
					type="response",
					single.tree=FALSE)
}

if(response.type=="categorical"){
	vote<-pred[,,1]
	pred<-colnames(vote)[apply(vote,1,which.max)]
	SGB.PRED<-data.frame(qdata.y,pred,vote)
	names(SGB.PRED)<-c("obs","pred",colnames(vote))
}else{
	SGB.PRED<-data.frame(obs=qdata.y,pred=pred)
}

rownames(SGB.PRED)<-rownames(qdata)

return(SGB.PRED)
}


#############################################################################################
########################## RF - Model Creation - Binary response ############################
#############################################################################################

rF.binary<-function(	qdata,
				predList,
				response.name,
				ntree=500,
				mtry=NULL,
				replace=TRUE,
				strata=NULL,
				sampsize = NULL,
				proximity = NULL,
				seed=NULL){

## This function generates a presence/absence (binary categorical) model using Random Forests.
##	Inputs: Full dataset, training indices, predictor names, response name, and seed (optional)
##	Output: Random Forest model


if(!is.null(seed)){
	set.seed(seed)}

qdata.x<-qdata[,match(predList,names(qdata))]

is.fact<-sapply(qdata.x,is.factor)
if(any(is.fact)){
	qdata.x[,is.fact]<-lapply(qdata.x[,is.fact,drop=FALSE],factor)}

qdata.y<-qdata[,response.name]
if(!is.numeric(qdata.y)){
	stop("If 'response.type is 'Binary' then 'response.name' must be numeric")}
qdata.y[qdata.y>0]<-1
qdata.y[qdata.y<0]<-0
qdata.y<-as.factor(qdata.y)

#print("about to start tuning")

if(is.null(mtry)){

	A<-list(	x=quote(qdata.x), y=quote(qdata.y),
			doBest=FALSE,
			importance=TRUE,
			proximity=FALSE,
			plot=FALSE,
			replace=replace,
			strata=strata,
			sampsize=sampsize)

	A<-A[!sapply(A, is.null)]

	RT<-do.call("tuneRF", A)

	mtry<-RT[which.min(RT[,2]),1]
}

#print("finished tuning")

A<-list(	x=quote(qdata.x), y=quote(qdata.y),
		importance=TRUE,
		proximity=proximity,
		mtry=mtry,
		ntree=ntree,
		replace=replace,
		strata=strata,
		sampsize=sampsize)

A<-A[!sapply(A, is.null)]

RF<-do.call("randomForest", A)

#rast.factors<-names(is.fact[is.fact])

RF$response<-response.name

return(RF)
}



#############################################################################################
############################ RF - Predict - Binary response #################################
#############################################################################################

prediction.rF.binary<-function(	prediction.type,
						qdata,
						response.name=RF$response,
						RF
						){

## This function makes predictions to test data for Random Forest presence/absence (binary 
## categorical) model.
##	Inputs: Training data, training data indices, predictor names, response variable name,
##			the Random Forest model and what to do if NAs are in predictors (default).
##	Output: Observed and predicted values.

if(is.null(response.name)){
	stop("must provide response name")}

if(prediction.type=="OOB"){
	pred<-predict(RF, type="vote")[,"1"]
	qdata.y<-RF$y
}
if(prediction.type=="TEST"){
	predList<-row.names(RF$importance)

	qdata.x<-qdata[,match(predList,names(qdata))]

	qdata.y<-qdata[,response.name]
	if(!is.numeric(qdata.y)){
		stop("If 'response.type is 'Binary' then 'response.name' must be numeric")}
	qdata.y[qdata.y>0]<-1
	qdata.y[qdata.y<0]<-0
	qdata.y<-as.factor(qdata.y)

	pred<-predict(RF, qdata.x,type="vote")[,"1"]
}

RF.PRED<-data.frame(	cbind(obs=as.numeric(as.character(qdata.y)), pred=pred))
rownames(RF.PRED)<-rownames(qdata)

return(RF.PRED)
}

#############################################################################################
######################## RF - Model Creation - Categorical response #########################
#############################################################################################

rF.categorical<-function(	qdata,
					predList,
					response.name,
					ntree=500,
					mtry=NULL,
					replace=TRUE,
					strata=NULL,
					sampsize = NULL,
					proximity = NULL,
					seed=NULL){

## This function generates a presence/absence (binary categorical) model using Random Forests.
##	Inputs: Full dataset, training indices, predictor names, response name, and seed (optional)
##	Output: Random Forest model


if(!is.null(seed)){
	set.seed(seed)}

qdata.x<-qdata[,match(predList,names(qdata))]

is.fact<-sapply(qdata.x,is.factor)
if(any(is.fact)){
	qdata.x[,is.fact]<-lapply(qdata.x[,is.fact,drop=FALSE],factor)}

qdata.y<-qdata[,response.name]
if(!is.factor(qdata.y)){
		qdata.y<-as.factor(qdata.y)}

#print("about to start tuning")

if(is.null(mtry)){

	A<-list(	x=quote(qdata.x), y=quote(qdata.y),
			doBest=FALSE,
			importance=TRUE,
			proximity=FALSE,
			plot=FALSE,
			replace=replace,
			strata=strata,
			sampsize=sampsize)

	A<-A[!sapply(A, is.null)]

	RT<-do.call("tuneRF", A)

	mtry<-RT[which.min(RT[,2]),1]
}

#print("finished tuning")

A<-list(	x=quote(qdata.x), y=quote(qdata.y),
		importance=TRUE,
		proximity=proximity,
		mtry=mtry,
		ntree=ntree,
		replace=replace,
		strata=strata,
		sampsize=sampsize)

A<-A[!sapply(A, is.null)]

RF<-do.call("randomForest", A)

#rast.factors<-names(is.fact[is.fact])

RF$response<-response.name

return(RF)
}



#############################################################################################
####################### RF - Predict - Categorical response #################################
#############################################################################################

prediction.rF.categorical<-function(	prediction.type,
							qdata,
							response.name=RF$response,
							RF
							){

## This function makes predictions to test data for Random Forest presence/absence (binary 
## categorical) model.
##	Inputs: Training data, training data indices, predictor names, response variable name,
##			the Random Forest model and what to do if NAs are in predictors (default).
##	Output: Observed and predicted values.

if(is.null(response.name)){
	stop("must provide response name")}


if(prediction.type=="OOB"){
	pred<-predict(RF)
	vote<-predict(RF, type="vote")
	qdata.y<-RF$y
}
if(prediction.type=="TEST"){

	predList<-row.names(RF$importance)

	qdata.x<-qdata[,match(predList,names(qdata))]

	qdata.y<-qdata[,response.name]
	if(!is.factor(qdata.y)){
			qdata.y<-as.factor(qdata.y)}

	pred<-predict(RF, qdata.x)
	vote<-predict(RF, qdata.x,type="vote")

}

RF.PRED<-data.frame(	obs=qdata.y, pred=pred, vote)
names(RF.PRED)<-c("obs","pred",colnames(vote))

rownames(RF.PRED)<-rownames(qdata)

return(RF.PRED)
}

#############################################################################################
#################### RF - Model Creation - Continuous Response ##############################
#############################################################################################

rF.continuous<-function(	qdata,
					predList,
					response.name,
					ntree=500,
					mtry=NULL,
					replace=TRUE,
					strata=NULL,
					sampsize = NULL,
					proximity = NULL,
					seed=NULL){


## This function generates a continuous response model using Random Forests.
##	Inputs: Full dataset, training indices, predictor names, response name, and seed (optional)
##	Output: Random Forest model

if(!is.null(seed)){
	set.seed(seed)}

qdata.x<-qdata[,match(predList,names(qdata))]

is.fact<-sapply(qdata.x,is.factor)
if(any(is.fact)){
	qdata.x[,is.fact]<-lapply(qdata.x[,is.fact,drop=FALSE],factor)}

qdata.y<-qdata[,response.name]



if(is.null(mtry)){

	A<-list(	x=quote(qdata.x), y=quote(qdata.y),
			doBest=FALSE,
			importance=TRUE,
			proximity=FALSE,
			plot=FALSE,
			replace=replace,
			strata=strata,
			sampsize=sampsize)

	A<-A[!sapply(A, is.null)]

	RT<-do.call("tuneRF", A)

	mtry<-RT[which.min(RT[,2]),1]
}

A<-list(	x=quote(qdata.x), y=quote(qdata.y),
		importance=TRUE,
		proximity=proximity,
		mtry=mtry,
		ntree=ntree,
		replace=replace,
		strata=strata,
		sampsize=sampsize)

A<-A[!sapply(A, is.null)]

RF<-do.call("randomForest", A)

#is.fact<-sapply(qdata.x,is.factor)
#rast.factors<-names(is.fact[is.fact])

RF$response<-response.name

return(RF)

}


#############################################################################################
########################## RF - Predict - Continuous Response ###############################
#############################################################################################

prediction.rF.continuous<-function(	prediction.type,
						qdata,
						response.name=RF$response,
						RF
						){

## This function makes predictions to test data for Random Forest continuous model.
##	Inputs: Training data, training data indices, predictor names, response variable name,
##			the Random Forest model and what to do if NAs are in predictors (default).
##	Output: Observed and predicted values.


if(is.null(response.name)){
	stop("must provide response name")}

if(prediction.type=="OOB"){
	pred<-RF$predicted
	qdata.y<-RF$y
}
if(prediction.type=="TEST"){
	predList<-row.names(RF$importance)
	qdata.x<-qdata[,match(predList,names(qdata))]
	pred<-predict(RF, qdata.x)
	qdata.y<-qdata[,response.name]
}

RF.PRED<-data.frame(	cbind(obs=as.numeric(as.character(qdata.y)),pred=pred))
rownames(RF.PRED)<-rownames(qdata)

return(RF.PRED)
}

#############################################################################################
#################### QRF - Model Creation - Continuous Response ######################
#############################################################################################

rF.quantreg<-function(	qdata,
				predList,
				response.name,
				ntree=1000,
				mtry=NULL,
				#replace=TRUE,
				#strata=NULL,
				#sampsize = NULL,
				#proximity = NULL,
				seed=NULL,
				importance=FALSE,
				quantiles=quantiles){


## This function generates a continuous response model using Random Forests.
##	Inputs: Full dataset, training indices, predictor names, response name, and seed (optional)
##	Output: quantile regression Random Forest model

if(!is.null(seed)){
	set.seed(seed)}

qdata.x<-qdata[,match(predList,names(qdata))]

is.fact<-sapply(qdata.x,is.factor)
if(any(is.fact)){
	qdata.x[,is.fact]<-lapply(qdata.x[,is.fact,drop=FALSE],factor)}

qdata.y<-qdata[,response.name]

if(is.null(mtry)){
	mtry <- ceiling(ncol(qdata.x)/3)
}

#print(paste("mtry =",mtry))

A<-list(	x=quote(qdata.x), y=quote(qdata.y),
		#importance=TRUE,
		#proximity=proximity,
		mtry=mtry,
		ntree=ntree,
		#replace=replace,
		#strata=strata,
		#sampsize=sampsize,
		importance=importance,
		quantiles=quantiles
	)

A<-A[!sapply(A, is.null)]

#f<-quantregForest::quantregForest
f<-getExportedValue(ns=asNamespace("quantregForest"), name="quantregForest")
QRF<-do.call(f, A)

#is.fact<-sapply(qdata.x,is.factor)
#rast.factors<-names(is.fact[is.fact])

QRF$response<-response.name

return(QRF)

}


#############################################################################################
########################## QRF - Predict - Continuous Response ##############################
#############################################################################################

prediction.rF.quantreg<-function(	prediction.type,
					qdata,
					response.name=QRF$response,
					QRF,
					quantiles=quantiles,
					all=FALSE
						){

## This function makes predictions to test data for Random Forest continuous model.
##	Inputs: Training data, training data indices, predictor names, response variable name,
##			the Random Forest model and what to do if NAs are in predictors (default).
##	Output: Observed and predicted values.


if(is.null(response.name)){
	stop("must provide response name")}

if(prediction.type=="OOB"){
	pred<-predict(QRF,quantiles=quantiles,all=all)
	qdata.y<-QRF$y
}
if(prediction.type=="TEST"){
	predList<-row.names(QRF$importance)
	qdata.x<-qdata[,match(predList,names(qdata))]
	pred<-predict(QRF, newdata=qdata.x,quantiles=quantiles,all=all)
	qdata.y<-qdata[,response.name]
}

quantlab<-sprintf("P%02d", 100*quantiles)

QRF.PRED<-data.frame(	cbind(obs=as.numeric(as.character(qdata.y)),pred))
names(QRF.PRED)   <-c("obs",quantlab)
rownames(QRF.PRED)<-rownames(qdata)

return(QRF.PRED)
}

#############################################################################################
########################## CF - Model Creation - Binary response ############################
#############################################################################################

CF.binary<-function(	qdata,
				predList,
				response.name,
				subset, 
				weights,
				controls,
				xtrafo, 
				ytrafo,
				scores,
				seed=NULL){

## This function generates a presence/absence (binary categorical) model using Random Forests.
##	Inputs: Full dataset, training indices, predictor names, response name, and seed (optional)
##	Output: Random Forest model

if(!is.null(seed)){
	set.seed(seed)}

qdata.x<-qdata[,match(predList,names(qdata))]

qdata.x<-CF.int2num(qdata.x)

is.fact<-sapply(qdata.x,is.factor)
if(any(is.fact)){
	qdata.x[,is.fact]<-lapply(qdata.x[,is.fact,drop=FALSE],factor)}

qdata.y<-qdata[,response.name]
if(!is.numeric(qdata.y)){
	stop("If 'response.type is 'Binary' then 'response.name' must be numeric")}
qdata.y[qdata.y>0]<-1
qdata.y[qdata.y<0]<-0
qdata.y<-as.factor(qdata.y)

qdata.formula<-cbind(qdata.y,qdata.x)
names(qdata.formula) <- c(response.name,names(qdata.x))

FORMULA <- as.formula(paste(response.name,paste(predList,collapse=" + "),sep=" ~ "))

A<-list(	formula = FORMULA, 
		data = qdata.formula, 
		subset = subset, 
		weights = weights, 
       	controls = controls,
        	xtrafo = xtrafo, 
		ytrafo = ytrafo, 
		scores = scores)

A<-A[!sapply(A, is.null)]

f<-getExportedValue(ns=asNamespace("party"),name="cforest")
CF<-do.call(f, A)

return(CF)
}



#############################################################################################
############################ CF - Predict - Binary response #################################
#############################################################################################

prediction.CF.binary<-function(	prediction.type,
						qdata,
						response.name,
						CF
						){

## This function makes predictions to test data for Random Forest presence/absence (binary 
## categorical) model.
##	Inputs: Training data, training data indices, predictor names, response variable name,
##			the Random Forest model and what to do if NAs are in predictors (default).
##	Output: Observed and predicted values.

if(is.null(response.name)){
	stop("must provide response name")}

qdata<-CF.int2num(qdata)

if(prediction.type=="OOB"){
	pred<-CF.list2df(predict(CF, OOB=TRUE, type="prob"))[,2]
	qdata.y<-CF@responses@variables[[1]]
}
if(prediction.type=="TEST"){
	#predList<-colnames(object@data@get("input"))
	#qdata.x<-qdata[,match(predList,names(qdata))]

	qdata.y<-qdata[,response.name]
	if(!is.numeric(qdata.y)){
		stop("If 'response.type is 'Binary' then 'response.name' must be numeric")}
	qdata.y[qdata.y>0]<-1
	qdata.y[qdata.y<0]<-0
	qdata.y<-as.factor(qdata.y)

	pred<-CF.list2df(predict(CF, newdata=qdata,OOB=FALSE,type="prob"))[,2]
}

CF.PRED<-data.frame(	cbind(obs=as.numeric(as.character(qdata.y)), pred=pred))
names(CF.PRED)<-c("obs","pred")
rownames(CF.PRED)<-rownames(qdata)

return(CF.PRED)
}

#############################################################################################
######################## CF - Model Creation - Categorical response #########################
#############################################################################################

CF.categorical<-function(	qdata,
					predList,
					response.name,
					subset, 
					weights,
					controls,
					xtrafo, 
					ytrafo,
					scores,
					seed=NULL){

## This function generates a presence/absence (binary categorical) model using Random Forests.
##	Inputs: Full dataset, training indices, predictor names, response name, and seed (optional)
##	Output: Random Forest model


if(!is.null(seed)){
	set.seed(seed)}

qdata.x<-qdata[,match(predList,names(qdata))]

qdata.x<-CF.int2num(qdata.x)

is.fact<-sapply(qdata.x,is.factor)
if(any(is.fact)){
	qdata.x[,is.fact]<-lapply(qdata.x[,is.fact,drop=FALSE],factor)}

qdata.y<-qdata[,response.name]
if(!is.factor(qdata.y)){
		qdata.y<-as.factor(qdata.y)}

qdata.formula<-cbind(qdata.y,qdata.x)
names(qdata.formula) <- c(response.name,names(qdata.x))

FORMULA <- as.formula(paste(response.name,paste(predList,collapse=" + "),sep=" ~ "))

A<-list(	formula = FORMULA, 
		data = qdata.formula, 
		subset = subset, 
		weights = weights, 
       	controls = controls,
        	xtrafo = xtrafo, 
		ytrafo = ytrafo, 
		scores = scores)

A<-A[!sapply(A, is.null)]

f<-getExportedValue(ns=asNamespace("party"),name="cforest")
CF<-do.call(f, A)

return(CF)
}



#############################################################################################
####################### CF - Predict - Categorical response #################################
#############################################################################################

prediction.CF.categorical<-function(	prediction.type,
							qdata,
							response.name,
							CF
							){

## This function makes predictions to test data for Random Forest presence/absence (binary 
## categorical) model.
##	Inputs: Training data, training data indices, predictor names, response variable name,
##			the Random Forest model and what to do if NAs are in predictors (default).
##	Output: Observed and predicted values.

if(is.null(response.name)){
	stop("must provide response name")}

qdata<-CF.int2num(qdata)

if(prediction.type=="OOB"){

	pred<-predict(CF, OOB=TRUE, type="response")
	vote<-CF.list2df(predict(CF, OOB=TRUE, type="prob"))
	colnames(vote)<-sapply(strsplit(colnames(vote),paste(response.name,".",sep="")),'[',2)

	qdata.y<-CF@responses@variables[[1]]
}
if(prediction.type=="TEST"){

	#predList<-colnames(object@data@get("input"))
	#qdata.x<-qdata[,match(predList,names(qdata))]

	qdata.y<-qdata[,response.name]
	if(!is.factor(qdata.y)){
			qdata.y<-as.factor(qdata.y)}

	pred<-predict(CF, newdata=qdata, OOB=FALSE, type="response")
	vote<-CF.list2df(predict(CF, newdata=qdata, OOB=FALSE, type="prob"))
	colnames(vote)<-sapply(strsplit(colnames(vote),paste(response.name,".",sep="")),'[',2)

}

CF.PRED<-data.frame(	obs=qdata.y, pred=pred, vote)
names(CF.PRED)<-c("obs","pred",colnames(vote))
rownames(CF.PRED)<-rownames(qdata)

return(CF.PRED)
}

#############################################################################################
#################### CF - Model Creation - Continuous Response ##############################
#############################################################################################

CF.continuous<-function(	qdata,
					predList,
					response.name,
					subset, 
					weights,
					controls,
					xtrafo, 
					ytrafo,
					scores,
					seed=NULL){


## This function generates a continuous response model using Random Forests.
##	Inputs: Full dataset, training indices, predictor names, response name, and seed (optional)
##	Output: Random Forest model

if(!is.null(seed)){
	set.seed(seed)}

qdata.x<-qdata[,match(predList,names(qdata))]

qdata.x<-CF.int2num(qdata.x)

is.fact<-sapply(qdata.x,is.factor)
if(any(is.fact)){
	qdata.x[,is.fact]<-lapply(qdata.x[,is.fact,drop=FALSE],factor)}

qdata.y<-qdata[,response.name]

qdata.formula<-cbind(qdata.y,qdata.x)
names(qdata.formula) <- c(response.name,names(qdata.x))

FORMULA <- as.formula(paste(response.name,paste(predList,collapse=" + "),sep=" ~ "))

A<-list(	formula = FORMULA, 
		data = qdata.formula, 
		subset = subset, 
		weights = weights, 
       	controls = controls,
        	xtrafo = xtrafo, 
		ytrafo = ytrafo, 
		scores = scores)

A<-A[!sapply(A, is.null)]

f<-getExportedValue(ns=asNamespace("party"),name="cforest")
CF<-do.call(f, A)

return(CF)
}


#############################################################################################
########################## CF - Predict - Continuous Response ###############################
#############################################################################################

prediction.CF.continuous<-function(	prediction.type,
						qdata,
						response.name,
						CF
						){

## This function makes predictions to test data for Random Forest continuous model.
##	Inputs: Training data, training data indices, predictor names, response variable name,
##			the Random Forest model and what to do if NAs are in predictors (default).
##	Output: Observed and predicted values.


if(is.null(response.name)){
	stop("must provide response name")}

qdata<-CF.int2num(qdata)

if(prediction.type=="OOB"){
	pred<-predict(CF, OOB=TRUE, type="response")
	qdata.y<-CF@responses@variables[[1]]
}
if(prediction.type=="TEST"){
	#predList<-colnames(object@data@get("input"))
	#qdata.x<-qdata[,match(predList,names(qdata))]

	pred<-predict(CF, newdata=qdata, OOB=FALSE, type="response")
	qdata.y<-qdata[,response.name]
}

CF.PRED<-data.frame(	cbind(obs=as.numeric(as.character(qdata.y)),pred=pred))
names(CF.PRED)<-c("obs","pred")
rownames(CF.PRED)<-rownames(qdata)

return(CF.PRED)
}


#############################################################################################
######################## Model Creation - Wrapper Function ##################################
#############################################################################################

create.model<-function(	qdata,
				model.type=NULL,		# "RF", "QRF", "SGB"
				folder=NULL,		# No ending slash, to output to working dir = getwd()
				predList,
				response.name=NULL,
				response.type,			# "binary", "continuous",
				seed=NULL,
				keep.data=TRUE,

			# RF arguments:
				ntree=500,
				mtry=NULL,
				replace=TRUE,
				strata=NULL,
				sampsize = NULL,
				proximity=proximity,

			# QRF arguments:
				importance=FALSE,
				quantiles=c(0.1,0.5,0.9),

			# CF arguments:
				subset=NULL, 
				weights=NULL,
				controls=party::cforest_unbiased(),
				xtrafo=party::ptrafo, 
				ytrafo=party::ptrafo,
				scores=NULL,

			# SGB arguments:
				n.trees=NULL,                 # number of trees
				shrinkage=0.001,   	      # shrinkage or learning rate,
                  	interaction.depth=10,		# 1: additive model, 2: two-way interactions, etc.
				bag.fraction = 0.5,          	# subsampling fraction, 0.5 is probably best
				nTrain = NULL,
				#train.fraction = NULL,       # fraction of data for training,
                  	n.minobsinnode = 10,         	# minimum total weight needed in each node
				var.monotone = NULL

){


### Set Seed ###
if(!is.null(seed)){
	set.seed(seed)}

if(model.type=="RF"){
	if(response.type=="binary"){
		#print("calling rF.binary")
		model.obj<-rF.binary(	qdata=qdata,
						predList=predList,
						response.name=response.name,
						ntree=ntree,
						mtry=mtry,
						replace=replace,
						strata=strata,
						sampsize=sampsize,
						proximity=proximity,
						seed=NULL)}
	if(response.type=="categorical"){
		#print("calling rF.categorical")
		model.obj<-rF.categorical(	qdata=qdata,
							predList=predList,
							response.name=response.name,
							ntree=ntree,
							mtry=mtry,
							replace=replace,
							strata=strata,
							sampsize=sampsize,
							proximity=proximity,
							seed=NULL)}
	if(response.type=="continuous"){
		#print("calling rF.continuous")
		model.obj<-rF.continuous(	qdata=qdata,
							predList=predList,
							response.name=response.name,
							ntree=ntree,
							mtry=mtry,
							replace=replace,
							strata=strata,
							sampsize=sampsize,
							proximity=proximity,
							seed=NULL)}
}

if(model.type=="QRF"){
		#print("creating quantreg model now: calling rF.quantreg")
		QRF<-rF.quantreg(	qdata=qdata,
						predList=predList,
						response.name=response.name,
						ntree=ntree,
						mtry=mtry,
						#replace=replace,
						#strata=strata,
						#sampsize=sampsize,
						#proximity=proximity,
						seed=NULL,
						importance=importance,
						quantiles=quantiles)
		#print("creating quantreg model now: calling rF.continuous")
		RF<-rF.continuous(	qdata=qdata,
					predList=predList,
					response.name=response.name,
					ntree=ntree,
					mtry=mtry,
					#replace=replace,
					#strata=strata,
					#sampsize=sampsize,
					proximity=proximity,
					seed=NULL)
		model.obj<-list(QRF=QRF,RF=RF)

}

if(model.type=="CF"){
	if(response.type=="binary"){
		#print("calling CF.binary")
		model.obj<-CF.binary(	qdata=qdata,
						predList=predList,
						response.name=response.name,
						subset=subset, 
						weights=weights,
						controls=controls,
						xtrafo=xtrafo, 
						ytrafo=ytrafo,
						scores=scores,
						seed=NULL)}
	if(response.type=="categorical"){
		#print("calling CF.categorical")
		model.obj<-CF.categorical(	qdata=qdata,
							predList=predList,
							response.name=response.name,
							subset=subset, 
							weights=weights,
							controls=controls,
							xtrafo=xtrafo, 
							ytrafo=ytrafo,
							scores=scores,
							seed=NULL)}
	if(response.type=="continuous"){
		#print("calling CF.continuous")
		model.obj<-CF.continuous(	qdata=qdata,
							predList=predList,
							response.name=response.name,
							subset=subset, 
							weights=weights,
							controls=controls,
							xtrafo=xtrafo, 
							ytrafo=ytrafo,
							scores=scores,
							seed=NULL)}
}

if(model.type=="SGB"){

	#print("calling model.SGB")
	model.obj<-model.SGB(	qdata=qdata,
					predList=predList,
					response.name=response.name,
					seed=NULL,
					response.type=response.type, 				
					n.trees=n.trees,                 	# number of trees
					shrinkage=shrinkage,   	      	# shrinkage or learning rate,
                  		interaction.depth=interaction.depth,# 1: additive model, 2: two-way interactions, etc.
					bag.fraction=bag.fraction,          # subsampling fraction, 0.5 is probably best
					nTrain=nTrain,
					#train.fraction=train.fraction,      # fraction of data for training,
                  		n.minobsinnode=n.minobsinnode,      # minimum total weight needed in each node
					keep.data=keep.data,
					var.monotone = var.monotone)

}

return(model.obj)
}

#############################################################################################
######################## Model Prediction - Wrapper Function ##################################
#############################################################################################




prediction.model<-function(	model.obj,
					model.type,
					qdata,
					folder=NULL,		# No ending slash, to output to working dir = getwd()
					response.name,
					response.type,
					unique.rowname="row_index",	# Row identifier

				# Model Evaluation Arguments
					prediction.type=NULL,
					MODELpredfn=NULL,
					v.fold=FALSE,

					na.action=na.action,			
					NA.ACTION=NA.ACTION,
					model.na.action=model.na.action,			
					model.NA.ACTION=model.NA.ACTION,

					model.levels=NULL,

				# QRF arguments
					quantiles = c(0.1, 0.5, 0.9),
					all=all,
					LOWER=0.10,
					UPPER=0.90,

	
				# CF arguments (prediction.type CV only)
					subset, 
					weights,
					controls,
					xtrafo, 
					ytrafo,
					scores,

				# SGB arguments
					n.trees
){

#should warnings be immeadiate
if(model.type=="CF" && prediction.type=="CV"){
	WARN.IM<-TRUE
}else{
	WARN.IM<-FALSE
}

#print("prediction wrapper function")
### make predictions ###
if(prediction.type!="CV"){
	if(model.type=="RF"){
			
		if(response.type=="binary"){
			#print("starting binary predictions")
			PRED<-prediction.rF.binary(	prediction.type=prediction.type,
								qdata=qdata,
								response.name=response.name,
								RF=model.obj)}
		if(response.type=="categorical"){
			#print("starting cat categorical predictions")
			PRED<-prediction.rF.categorical(	prediction.type=prediction.type,
									qdata=qdata,
									response.name=response.name,
									RF=model.obj)}
		if(response.type=="continuous"){
			#print("starting continuous predictions")
			PRED<-prediction.rF.continuous(	prediction.type=prediction.type,
									qdata=qdata,
									response.name=response.name,
									RF=model.obj)}
	}

	if(model.type=="QRF"){

		if("QRF"%in%names(model.obj)){
			#print("starting continuous predictions")
			PRED<-prediction.rF.quantreg(	prediction.type=prediction.type,
									qdata=qdata,
									response.name=response.name,
									QRF=model.obj$QRF,
									quantiles=c(quantiles,LOWER,0.50,UPPER),
									all=all)
			RFPRED<-prediction.rF.continuous(	prediction.type=prediction.type,
									qdata=qdata,
									response.name=response.name,
									RF=model.obj$RF)
			Nq<-length(quantiles)
			PRED<-cbind(PRED[,1:(Nq+1)],RFPRED$pred,PRED[,(Nq+2):(Nq+4)])
			names(PRED)<-c(names(PRED[,1:(Nq+1)]),"RFmean","lower","median","upper")
		}else{
			PRED<-prediction.rF.quantreg(	prediction.type=prediction.type,
									qdata=qdata,
									response.name=response.name,
									QRF=model.obj,
									quantiles=c(quantiles,LOWER,0.50,UPPER),
									all=all)
			Nq<-length(quantiles)
			names(PRED)<-c(names(PRED[,1:Nq+1]),"lower","median","upper")
		}
	}

	if(model.type=="CF"){	
		if(response.type=="binary"){
			#print("starting binary predictions")
			PRED<-prediction.CF.binary(	prediction.type=prediction.type,
								qdata=qdata,
								response.name=response.name,
								CF=model.obj)}
		if(response.type=="categorical"){
			#print("starting cat categorical predictions")
			PRED<-prediction.CF.categorical(	prediction.type=prediction.type,
									qdata=qdata,
									response.name=response.name,
									CF=model.obj)}
		if(response.type=="continuous"){
			#print("starting continuous predictions")
			PRED<-prediction.CF.continuous(	prediction.type=prediction.type,
									qdata=qdata,
									response.name=response.name,
									CF=model.obj)}
	}


	if(model.type=="SGB"){
		#print("calling prediction.SGB")
		PRED<-prediction.SGB(	prediction.type=prediction.type,
						qdata=qdata,
						response.name=response.name,
						SGB=model.obj,
						n.trees=n.trees)
	}

	#print("got predictions, checking rownames")

	#If prediction type is OOB, model was roughfix and diagnostics is omit, must remove na rows from predictions
	#if(!is.na(NA.ACTION) && !is.na(model.NA.ACTION)){
		if(prediction.type=="OOB" && model.NA.ACTION%in%"roughfix" && NA.ACTION%in%"omit"){

			#print(paste("     nrow(qdata) =",length(rownames(qdata))))
			#print(paste("          remove these rows =",paste(model.obj$na.action,collapse=",")))

			#print(paste("     nrow(PRED) =", nrow(PRED)))
			PRED<-PRED[-model.obj$na.action,] 
			#print(paste("     nrow(PRED) =", nrow(PRED)))
		}
	#}


	#create qdata.y that has been transformed as response gets transformed

	qdata.y<-qdata[,response.name]
	
	if(response.type=="binary"){
		if(!is.numeric(qdata.y)){
			stop("If 'response.type is 'Binary' then 'response.name' must be numeric")}
		qdata.y[qdata.y>0]<-1
		qdata.y[qdata.y<0]<-0
	}
	if(response.type=="categorical"){
		if(!is.factor(qdata.y)){qdata.y<-as.factor(qdata.y)}
	}
	#print(paste("levels qdata.y =",levels(qdata.y)))	

	#If prediction is OOB check that training data responses are the same as those used to build the model
	if(prediction.type=="OOB"){
		if(!isTRUE(all.equal(PRED$obs,qdata.y))){
			stop("Prediction type is 'OOB' but responses in 'qdata.trainfn' do not match data used to build the model")
		}
	}

	#Note - If this warning appears there is a bug
	if(!isTRUE(all.equal(PRED$obs,qdata.y))){
		stop("Something has gone wrong and observed responses are scrambled")
	}

	rownames(PRED)<-rownames(qdata)

}else{
	print(paste("Begining ",v.fold,"-fold cross validation:",sep=""))

	predList<-check.predList(model.obj=model.obj,model.type=model.type)

	if(model.type=="RF"){
		ntree<-model.obj$ntree
		mtry<-model.obj$mtry
		replace<-model.obj$call$replace}
	if(model.type=="QRF"){
		ntree  <-if("QRF"%in%names(model.obj)){model.obj$QRF$ntree}else{model.obj$ntree}
		mtry   <-if("QRF"%in%names(model.obj)){model.obj$QRF$mtry }else{model.obj$mtry}
		replace<-if( "RF"%in%names(model.obj)){model.obj$RF$call$replace}else{TRUE}
	}
	#if(model.type=="CF"){
	#	weights<-slotNames(model.obj@weights)
	#	controls
	#	xtrafo
	#	ytrafo
	#	scores
	#	}
	if(model.type=="SGB"){
		shrinkage<-model.obj$shrinkage
		interaction.depth<-model.obj$interaction.depth
		bag.fraction<-model.obj$bag.fraction
		nTrain<-model.obj$nTrain
		n.minobsinnode<-model.obj$n.minobsinnode
		
		#deal with nTrain<nrow(qdata)
		#print(paste("nTrain =",nTrain))
		#print(paste("nrow(qdata) =",nrow(qdata)))
		if(nTrain<nrow(qdata)){
			qdata<-qdata[1:nTrain,]}	
	}

	n.data=nrow(qdata)
	n.per.fold<-floor(n.data/v.fold)
	cv.index<-sample(rep(1:v.fold,(n.per.fold+1))[1:n.data])

	print(table(cv.index,qdata[,response.name]))
	
	###for binary models, make sure response is 0/1 that has been transformed as response gets transformed###

	if(response.type=="binary"){
		if(!is.numeric(qdata[,response.name])){
			stop("If 'response.type is 'Binary' then 'response.name' must be numeric")}
		qdata[,response.name][qdata[,response.name]>0]<-1
		qdata[,response.name][qdata[,response.name]<0]<-0
	}


	###cv.index[qdata[,p]=="70"]<-2###
	
	if(model.type=="RF" || model.type=="CF" || model.type=="SGB"){
		if(response.type=="categorical"){
			qdata.y<-qdata[,response.name]
			if(!is.factor(qdata.y)){qdata.y<-as.factor(qdata.y)}
			#print(paste("levels qdata.y =",levels(qdata.y)))
			PRED<-data.frame(matrix(0,0,3+length(levels(qdata.y))))
			names(PRED)<-c("obs","pred",levels(qdata.y),"Vfold")
		}else{
			PRED<-data.frame(matrix(0,0,2))
			names(PRED)<-c("obs","pred")
		}
	}
	if(model.type=="QRF"){
		PRED<-data.frame(matrix(0,0,5+length(quantiles)))
		names(PRED)<-c("obs",paste("quantile=", quantiles,sep=""),"RFmean","lower","median","upper")
	}


	###start folds###
	for(i in 1:v.fold){
		print(paste("starting fold", i))
		train.cv<-(1:nrow(qdata))[cv.index!=i]
		qdata.train.cv<-qdata[train.cv,]
		qdata.test.cv<-qdata[-train.cv,]

		###check for factor levels not found in other folds###
		print("starting factor checks")

		#print("Response Levels - Data")
		#print(levels(qdata.y))
		#print("Response levels - CV")
		#print(levels(qdata.train.cv[,response.name]))

		if(!is.null(model.levels)){
			predFactor<-names(model.levels)
			invalid.levels<-vector("list", length=length(predFactor))
			names(invalid.levels)<-predFactor

			for(p in predFactor){

				#train.levels<-levels(qdata.train.cv[,p])[levels(qdata.train.cv[,p])%in%qdata.train.cv[,p]]
				#test.levels<-levels(qdata.test.cv[,p])[levels(qdata.test.cv[,p])%in%qdata.test.cv[,p]]

				qdata.train.cv[,p]<-as.character(qdata.train.cv[,p])
				qdata.test.cv[,p]<-as.character(qdata.test.cv[,p])

				train.levels<-sort(unique(qdata.train.cv[,p]))
				test.levels<- sort(unique(qdata.test.cv[,p] ))

				invalid<-test.levels[!(test.levels%in%c(train.levels,NA))] 
				invalid.levels[[p]]<-invalid

				qdata.train.cv[,p]<-factor(qdata.train.cv[,p],levels=train.levels)
				qdata.test.cv[,p] <-factor(qdata.test.cv[,p],levels=train.levels)

				print(paste("     predictor            =",p))
				print(paste("     train levels[p]      =",paste(train.levels,collapse=",")))
				print(paste("     test levels[p]       =",paste(test.levels,collapse=",")))
				#print(paste("     invalid              =",paste(invalid,collapse=",")))
				#print(paste("     invalid levels[[p]]  =",paste(invalid.levels[[p]],collapse=",")))
				#print(paste("          invalid length =",length(invalid.levels[[p]])))
				#print(paste("     levels qdata.test.cv =",paste(levels(qdata.test.cv[,p]),collapse=",")))

				if(length(invalid.levels[[p]])>0){
					#print("     NA action triggered")
					N.invalid<-sum(is.na(qdata.test.cv[,p]))
					warn.lev<-paste(invalid.levels[[p]],collapse=", ")
					warning(	"Factored predictor ",p," in fold ",i, 
							" contains ", N.invalid, " data points of levels ",warn.lev,
							" not found in other folds and these levels treated as NA",
							immediate. = WARN.IM)
				}
			}
			
			NA.pred<-apply(qdata.test.cv[,predList],1,function(x){any(is.na(x))})
			#print(paste(     "NA.pred =",any(NA.pred)))
			#print("TEST TEST TEST")
			if(any(NA.pred)){
				#print("about to treat na's")
				#print("NA.ACTION:")
				#print(NA.ACTION)
				#print("na.action:")
				#print(na.action)
				qdata.test.cv<-na.action(qdata.test.cv)
				#print("finished treating NA's")
			}
			#print(qdata.test.cv[,predList])
		}

		###build models and make predictions###
		if(model.type=="RF"){
			if(response.type=="binary"){
				RF.cv<-rF.binary(	qdata=qdata.train.cv,
							predList=predList,
							response.name=response.name,
							ntree=ntree,
							mtry=mtry,
							replace=replace,
							seed=NULL)
				
				#print(paste("        calling prediction.rF.binary for fold",i))
				PRED.cv<-prediction.rF.binary(	prediction.type="TEST",
										qdata=qdata.test.cv,
										response.name=response.name,
										RF=RF.cv)
			}
			if(response.type=="categorical"){
				RF.cv<-rF.categorical(	qdata=qdata.train.cv,
							predList=predList,
							response.name=response.name,
							ntree=ntree,
							mtry=mtry,
							replace=replace,
							seed=NULL)
				#print(paste("        calling prediction.rF.categorical for fold",i))
				PRED.cv<-prediction.rF.categorical(	prediction.type="TEST",
										qdata=qdata.test.cv,
										response.name=response.name,
										RF=RF.cv)
			}

			if(response.type=="continuous"){
				##("generating new model")
				RF.cv<-rF.continuous(	qdata=qdata.train.cv,
								predList=predList,
								response.name=response.name,
								ntree=ntree,
								mtry=mtry,
								replace=replace,
								seed=NULL)
				#print(paste("        calling prediction.rF.continuous for fold",i))
				PRED.cv<-prediction.rF.continuous(	prediction.type="TEST",
										qdata=qdata.test.cv,
										response.name=response.name,
										RF=RF.cv)
			}
		}

		if(model.type=="QRF"){

			#print(paste("qdata:        ",qdata.train.cv[1,]))
			#print(paste("predList:     ",predList))
			#print(paste("response.name:",response.name))
			#print(paste("ntree:        ",ntree))
			#print(paste("mtry:         ",mtry))
			#print(paste("seed:         ",seed))
			
			QRF.cv<-rF.quantreg(	qdata=qdata.train.cv,
							predList=predList,
							response.name=response.name,
							ntree=ntree,
							mtry=mtry,
							#replace=replace,
							seed=NULL,
							importance=FALSE,
							quantiles=NULL)
			#print(paste("        calling prediction.rF.quantreg for fold",i))
			QRFPRED.cv<-prediction.rF.quantreg( prediction.type="TEST",
									qdata=qdata.test.cv,
									response.name=response.name,														
									QRF=QRF.cv,
									quantiles=c(quantiles,0.05,0.50,0.95),
									all=all)

			RF.cv<-rF.continuous(	qdata=qdata.train.cv,
						predList=predList,
						response.name=response.name,
						ntree=ntree,
						mtry=mtry,
						replace=replace,
						seed=NULL)
			RFPRED.cv<-prediction.rF.continuous(prediction.type="TEST",
								qdata=qdata.test.cv,
								response.name=response.name,
								RF=RF.cv)
			Nq<-length(quantiles)
			PRED.cv<-cbind(QRFPRED.cv[,1:(Nq+1)],RFPRED.cv$pred,QRFPRED.cv[,(Nq+2):(Nq+4)])
			names(PRED.cv)<-c(names(QRFPRED.cv[,1:(Nq+1)]),"RFmean","lower","median","upper")
		}

		if(model.type=="CF"){
			if(response.type=="binary"){
				CF.cv<-CF.binary(	qdata=qdata.train.cv,
							predList=predList,
							response.name=response.name,
							subset=subset, 
							weights=weights,
							controls=controls,
							xtrafo=xtrafo, 
							ytrafo=ytrafo,
							scores=scores,
							seed=NULL)
				
				
				#print(paste("        calling prediction.CF.binary for fold",i))
				PRED.cv<-prediction.CF.binary(	prediction.type="TEST",
										qdata=qdata.test.cv,
										response.name=response.name,
										CF=CF.cv)
			}
			if(response.type=="categorical"){
				CF.cv<-CF.categorical(	qdata=qdata.train.cv,
								predList=predList,
								response.name=response.name,
								subset=subset, 
								weights=weights,
								controls=controls,
								xtrafo=xtrafo, 
								ytrafo=ytrafo,
								scores=scores,
								seed=NULL)
				#print(paste("        calling prediction.CF.categorical for fold",i))
				PRED.cv<-prediction.CF.categorical(	prediction.type="TEST",
										qdata=qdata.test.cv,
										response.name=response.name,
										CF=CF.cv)
			}

			if(response.type=="continuous"){
				##("generating new model")
				CF.cv<-CF.continuous(	qdata=qdata.train.cv,
								predList=predList,
								response.name=response.name,
								subset=subset, 
								weights=weights,
								controls=controls,
								xtrafo=xtrafo, 
								ytrafo=ytrafo,
								scores=scores,
								seed=NULL)
				#print(paste("        calling prediction.CF.continuous for fold",i))
				PRED.cv<-prediction.CF.continuous(	prediction.type="TEST",
										qdata=qdata.test.cv,
										response.name=response.name,
										CF=CF.cv)
			}
		}

		if(model.type=="SGB"){
			#print(paste("      making SGB model for fold",i))

			
			SGB.cv<-suppressWarnings(model.SGB(	qdata=qdata.train.cv,
							predList=predList,
							response.name=response.name,
							seed=NULL,
							response.type=response.type, 				
							n.trees=n.trees,	               	# number of trees
							shrinkage=shrinkage,   	      	# shrinkage or learning rate,
                  				interaction.depth=interaction.depth,# 1: additive model, 2: two-way interactions, etc.
							bag.fraction=bag.fraction,          # subsampling fraction, 0.5 is probably best
							nTrain=nrow(qdata.train.cv),
							#train.fraction=train.fraction,      # fraction of data for training,
                  				n.minobsinnode=n.minobsinnode       # minimum total weight needed in each node
							))
			#print(paste("      making SGB predictions for fold",i))
			PRED.cv<-prediction.SGB(	prediction.type="TEST",
								qdata=qdata.test.cv,
								response.name=response.name,
								SGB=SGB.cv,
								n.trees=n.trees)
		}
		PRED.cv$VFold<-i
		PRED<-rbind(PRED,PRED.cv)
		#print(paste("     ending fold",i))
	}
PRED<-PRED[match(row.names(qdata),row.names(PRED),nomatch=0),] #cv only
}

PRED<-cbind(rownames(PRED),PRED)
colnames(PRED)[1]<-unique.rowname

PREDICTIONfn<-paste(MODELpredfn,".csv",sep="")
write.table(PRED,file=PREDICTIONfn,sep=",",row.names=FALSE)

return(PRED)
}

#############################################################################################
#############################################################################################
###################### Production Prediction - Sub Functions ################################
#############################################################################################
#############################################################################################



#############################################################################################
################################## Check filename ###########################################
#############################################################################################


FNcheck<-function(OUTPUTfn,folder,ERROR.NAME){
### checks filename for valid extensions, stops if there is more than 1 '.' in filename or invalid extension after the dot.
### If no, extension, warns and adds default extension of .img
### If no path, adds path from 'folder'.

	validExtensions<-c(".tif",
				".tiff",
				".grd",
				".asc",
				".nc",
				".cdf",
				".ncdf",
				".kml",
				".kmz",
				".big",
				".sgrd",
				".sdat",
				".bil",
				".bsq",
				".bip",
				".bmp",
				".gen",
				".bt",
				".envi",
				".ers",
				".img",
				".rst",
				".mpr",
				".rsw") 


	OUTPUTsplit<-strsplit(basename(OUTPUTfn),split="\\.")[[1]]
	OUTPUText<-paste(".",tail(OUTPUTsplit, 1),sep="")

	#if basename of output filename has more than 1 '.'
	if(length(OUTPUTsplit) > 2){
		stop(ERROR.NAME," ",OUTPUTfn," has more than one '.' The ",ERROR.NAME," should only use the character '.' to indicate the file extension.")}

	#if basename has exactly one dot, check if extension is valid
	if(length(OUTPUTsplit) == 2){
		if(!tolower(OUTPUText)%in%validExtensions){
			#OUTPUTfn<-paste(dirname(OUTPUTfn),"/",OUTPUTsplit[1],".img",sep="")
			stop(ERROR.NAME," extension ",OUTPUText," is invalid. See 'writeFormats()' for a list of valid raster file types.")}}

	#if basename has no extension, add default extension of '.img'
	if(length(OUTPUTsplit) == 1){
		OUTPUTfn<-paste(OUTPUTfn,".img",sep="")
		warning(ERROR.NAME," ",OUTPUTsplit," does not include an extension, using default extension of '.img'. New ",ERROR.NAME," is ",OUTPUTfn)
	}

	# if OUTPUTfn has no directory path, add folder to name
	if(identical(basename(OUTPUTfn),OUTPUTfn))
		{OUTPUTfn<-file.path(folder,OUTPUTfn)}

	return(OUTPUTfn)
}


#############################################################################################
################################# Fix projections ###########################################
#############################################################################################


projfix<-function(RAST,OUTPUTfn.noext){
#RAST	list of rasters suitable to be given to stack() function

	#find all projections
	PROJ<-lapply(RAST,projection)

	#create output file of all projections
	OUTPUTfn.proj<-paste(OUTPUTfn.noext,"_projections.txt",sep="")
	PROJout<-data.frame(predictor=names(PROJ),projection=sapply(PROJ,I))
	write.table(PROJout,file=OUTPUTfn.proj,row.names=FALSE,sep="\t")

	#find source filenames
	FN<-sapply(RAST,function(x){basename(filename(x))})

	#check if all rasters have projections
	HAVEPROJ<-sapply(RAST,function(x){is.na(projection(x)) || is.null(projection(x))})
	if(any(HAVEPROJ)){
		if(sum(HAVEPROJ)==1){
			warning.message<-paste("raster for predictor ",names(HAVEPROJ[HAVEPROJ])," is missing a projection, see '", OUTPUTfn.proj, "' for details", sep="")
		}else{
			warning.message<-paste("rasters for predictors ",paste(names(HAVEPROJ[HAVEPROJ]),collapse=" ")," are missing projections, see '", OUTPUTfn.proj, "' for details", sep="")
		}
		stop(warning.message)
	}


#	#check if projections are similar enough to reconcile
#	SAME<-sapply(PROJ,function(x,proj){projcompare(x,proj)},proj=PROJ[[1]])
#	SAME.R<-sapply(RAST,function(x,control.rast){compareRaster(x,control.rast,stopiffalse=FALSE)},control.rast=RAST[[1]])
#
#	if(all(SAME)){
#		if(all(SAME.R)){
#			#they all match, no need to do anything
#			RAST.out<-RAST
#		}else{
#			#set all projections to match layer 1
#			RAST.out<-RAST
#			RAST.out<-lapply(RAST,function(x,rast){projection(x)<-projection(rast); x},rast=RAST[[1]])
#			warning("one or more predictor layers had projections that needed to be reconciled, see '", OUTPUTfn.proj, "' for details")
#		}
#	}else{
#		stop("projections of predictor layers too different to reconcile, see '", OUTPUTfn.proj, "' for details")
#	}

	#check if projections are the same
	#SAME<-sapply(PROJ,function(x,proj){projcompare(x,proj)},proj=PROJ[[1]])
	SAME.R<-sapply(RAST,function(x,control.rast){compareRaster(x,control.rast,stopiffalse=FALSE)},control.rast=RAST[[1]])

	#if(all(SAME)){
	if(all(SAME.R)){
			#they all match, no need to do anything
			print("all predictor layer rasters match")
			RAST.out<-RAST
		#}else{
		#	#set all projections to match layer 1
		#	RAST.out<-RAST
		#	RAST.out<-lapply(RAST,function(x,rast){projection(x)<-projection(rast); x},rast=RAST[[1]])
		#	warning("one or more predictor layers had projections that needed to be reconciled, see '", OUTPUTfn.proj, "' for details")
		#}
	}else{
		stop("predictor layer rasters too different to reconcile, see '", OUTPUTfn.proj, "' for details")
	}

	return(RAST.out)
}

#############################################################################################
###################################### grd to gri ###########################################
#############################################################################################

grd2gri<-function(x){
	paste(strsplit(x,".grd")[[1]],".gri",sep="")}

#############################################################################################
###################################### Compare projections ##################################
#############################################################################################
# Note this function is not currently in use, it is currently replaced by raster package function compareRaster()

	projcompare <- function(layer1prj, layer2prj){
		########################################################################################
		## DESCRIPTION: Internal function to compare projections. If not the same, return FALSE
		##
		## ARGUMENTS:
		##   layer1prj - the projection name for layer1
		##   layer2prj - the projection name for layer1	
		##
		## NOTE:
		## Use function proj4string(layer) or projection(layer) to get the projection name 
		########################################################################################


		ellpsna <- FALSE
		datumna <- FALSE

		if(is.null(layer1prj) | is.null(layer2prj)){
			stop("check projection layers")
		}
	
		projnm1 <- strsplit(strsplit(layer1prj, "+proj=")[[1]][2], " +")[[1]][1]
		ellpsnm1 <- strsplit(strsplit(layer1prj, "+ellps=")[[1]][2], " +")[[1]][1]
		datumnm1 <- strsplit(strsplit(layer1prj, "+datum=")[[1]][2], " +")[[1]][1]
	
		projnm2 <- strsplit(strsplit(layer2prj, "+proj=")[[1]][2], " +")[[1]][1]
		ellpsnm2 <- strsplit(strsplit(layer2prj, "+ellps=")[[1]][2], " +")[[1]][1]
		datumnm2 <- strsplit(strsplit(layer2prj, "+datum=")[[1]][2], " +")[[1]][1]

		if(is.na(ellpsnm1) | is.na(ellpsnm2)){
			ellpsna <- TRUE
		}
		if(is.na(datumnm1) | is.na(datumnm2)){
			datumna <- TRUE
		}

		if(is.null(ellpsna) & is.null(datumna)){
			stop("CHECK PROJECTIONS. NEED ELLPSNM OR DATUM DEFINED IN BOTH PROJECTIONS")
		}else{
	
			## ASSUME WGS84 and NAD83 ARE EQUAL
			if(!is.na(datumnm1)){ if(datumnm1 == "WGS84"){ datumnm1 = "NAD83" } }
			if(!is.na(datumnm2)){ if(datumnm2 == "WGS84"){ datumnm2 = "NAD83" } }				

			if(projnm1 == projnm2){
				if(ellpsnm1 == ellpsnm2 | ellpsna){
					if(datumnm1 == datumnm2 | datumna){
						projmatch <- TRUE
					}else{
						projmatch <- FALSE
					}
				}else{
					projmatch <- FALSE
				}
			}else{
				projmatch <- FALSE
			}
		}
		return(projmatch)
	}

#############################################################################################
#############################################################################################
############################## color translation ############################################
#############################################################################################
#############################################################################################

col2trans<-function(col.names,alpha=0.5){
	col.out <- rgb(t(col2rgb(col.names)/255),alpha=alpha)
	return(col.out)
}

#############################################################################################
#############################################################################################
#################### Production Prediction - Actual Function ################################
#############################################################################################
#############################################################################################

production.prediction<-function(	model.obj,
						model.type,
						rastLUT,
						#na.action=NULL,
						NA.ACTION=NULL,
						NAval=-9999,
						model.levels,
						response.type,
						#response.name,
						keep.predictor.brick,							
						map.sd=FALSE,
						OUTPUTfn,
						OUTPUTfn.noext,
						#OUTPUTpath,
						OUTPUTname,
						OUTPUText,
						quantiles,
						n.trees){

#############################################################################################
################################## Create File names ########################################
#############################################################################################


### Creat filename for native raster format map output

TMPfn.map <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_map_",sep=""))

if(model.type=="QRF"){

	model.obj.QRF<-if("QRF"%in%names(model.obj)){model.obj$QRF}else{model.obj}
	model.obj.RF <-if( "RF"%in%names(model.obj)){model.obj$RF}else{NULL}

	OUTPUTfn  <- paste(OUTPUTfn.noext,"_QRF",OUTPUText,sep="")
	TMPfn.map <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_QRF_",sep=""))

	if(!is.null(model.obj.RF)){
		OUTPUTfn.RF<- paste(OUTPUTfn.noext,"_RF",OUTPUText,sep="")
		TMPfn.RF   <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_RF_",sep=""))
	}else{
		map.sd <- FALSE
	}
}


### Creat filename for predictor raster brick

if(keep.predictor.brick){
	OUTPUTfn.brick <- paste(OUTPUTfn.noext,"_brick",sep="")
}else{
	OUTPUTfn.brick <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_brick_",sep=""))
}

### map sd filenames

if(map.sd && model.type%in%c("RF","QRF") && response.type=="continuous"){

	OUTPUTfn.mean  <- paste(OUTPUTfn.noext,"_mean",OUTPUText,sep="")
	OUTPUTfn.stdev <- paste(OUTPUTfn.noext,"_stdev",OUTPUText,sep="")
	OUTPUTfn.coefv <- paste(OUTPUTfn.noext,"_coefv",OUTPUText,sep="")

	TMPfn.mean     <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_mean_",sep=""))
	TMPfn.stdev    <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_stdev_",sep=""))
	TMPfn.coefv    <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_coefv_",sep=""))	

}else{
	map.sd<-FALSE
}

#####################################################################################
########################## Extract predictor names ##################################
#####################################################################################


## Make sure raster predictor names match names in training/test data.

predList<-check.predList(model.obj=model.obj,model.type=model.type)

predLUT<-rastLUT[match(predList, rastLUT[,2]),] #only the predictors from the model, in same order as in model
if(any(predList!=predLUT[,2])){
	stop("predictor names from model do not match short names in rastLUT")}	

anyFactor<-FALSE
predFactor<-NULL
if(!is.null(model.levels)){
	anyFactor<-TRUE
	predFactor<-names(model.levels)
	extraLevels<-vector("list",length(predFactor))
	names(extraLevels)<-predFactor
}

print(paste("predFactor:",predFactor))

#####################################################################################
########################## Extract response classes #################################
#####################################################################################

if(response.type=="categorical"){
	if(model.type=="RF"){ Rclasses<-colnames(model.obj$votes)}
	#if(model.type=="CF"){ Rclasses<-model.obj@responses@levels[[response.name]]}#might be needed for multivariate?
	if(model.type=="CF"){ Rclasses<-model.obj@responses@levels[[1]]}
	if(model.type=="SGB"){Rclasses<-model.obj$classes}
	Nrc<-length(Rclasses)
}

########################################################################################
################################ Set Data Type #########################################
########################################################################################

if(response.type=="binary"){data.type <- "FLT4S"}
if(response.type=="categorical"){
	data.type <- "INT2S"
	Ylev<-Rclasses}
if(response.type=="continuous"){data.type <- "FLT4S"}

########################################################################################
########################### Build raster stack #########################################
########################################################################################

#require(raster)

RAST<-vector("list", 0)

for(p in predList){
	rastfn<-rastLUT[rastLUT[,2]==p,1]
	band<-  rastLUT[rastLUT[,2]==p,3]

	RAST[[p]]<-raster(rastfn,band=band)
}

RAST<-projfix(RAST,OUTPUTfn.noext=OUTPUTfn.noext)

#compareRaster(RAST,stopiffalse=TRUE, showwarning=TRUE) #, crs=FALSE)

RS<-stack(RAST)
RB<-brick(RS,values=TRUE,filename=OUTPUTfn.brick,overwrite=TRUE)

print("brick done")
########################################################################################
############################# Loop through rows ########################################
########################################################################################

print(paste("write start", OUTPUTfn))

#print(paste("datatype =",data.type))

if(model.type%in%c("RF","CF","SGB")){ 
	out <- raster(RAST[[1]])
	dataType(out) <- data.type
	NAvalue(out) <- NAval

	out <- writeStart(out, filename=TMPfn.map, overwrite=TRUE, datatype=data.type)
	#out <- writeStart(out, overwrite=TRUE,  dataType=data.type) #doesn't work
}
if(model.type == "QRF"){ 
	out <- brick(RAST[[1]],nl=length(quantiles),values=FALSE)
	names(out) <- paste("quantile=", quantiles)

	dataType(out) <- data.type
	NAvalue(out) <- NAval

	out <- writeStart(out, filename=TMPfn.map, overwrite=TRUE, datatype=data.type)
	#out <- writeStart(out, overwrite=TRUE,  dataType=data.type) #doesn't work

	if(!is.null(model.obj.RF)){
		out.RF <- raster(RAST[[1]])
		out.RF <- writeStart(out.RF, filename=extension(TMPfn.RF,""), overwrite=TRUE, datatype=data.type)}
}

if(map.sd){
	out.mean <- raster(RAST[[1]])
	out.mean <- writeStart(out.mean, filename=extension(TMPfn.mean,""), overwrite=TRUE, datatype=data.type)

	out.stdev <- raster(RAST[[1]])
	out.stdev <- writeStart(out.stdev, filename=extension(TMPfn.stdev,""), overwrite=TRUE, datatype=data.type)

	out.coefv <- raster(RAST[[1]])
	out.coefv <- writeStart(out.coefv, filename=extension(TMPfn.coefv,""), overwrite=TRUE, datatype=data.type)
}

###############################################################################################################
print("starting row by row predictions")
for(r in 1:(dim(RB)[1])){
	print(paste("     rows =",r))

	v <- data.frame(getValues(RB, r))

	v<-CF.int2num(v)

	### deal with factored predictors ###
	if(anyFactor){
		for(f in predFactor){
			f.lev<-model.levels[[f]]
			f.extra<-v[,f][!v[,f]%in%f.lev]
			extraLevels[[f]]<-unique(c(extraLevels[[f]],f.extra))
			
			#v[,f][!v[,f]%in%f.lev] <- ??? #leaving this line just in case we decide to add rough fix back in

			v[,f]<-factor(v[,f],levels=f.lev)
		}
	}

	### deal with NAval ###
	nonPredict <- apply(((v == NAval)|is.na(v)), 1, any)

	
	if(model.type%in%c("RF","CF","SGB")){
		v.pred<-rep(NA,length=nrow(v))
	}
	if(model.type=="QRF"){
		v.pred <- matrix(NA,nrow=nrow(v),ncol=length(quantiles))
		colnames(v.pred) <- paste("quantile=", quantiles)

		if(!is.null(model.obj.RF)){v.pred.RF<-rep(NA,length=nrow(v))}
	}
		
	if(any(!nonPredict)){
		if(model.type=="RF"){
			if(response.type=="binary"){
				v.pred[!nonPredict] <- signif(predict(model.obj, v[!nonPredict,,drop=FALSE],type="vote")[,"1"],2)}

			if(response.type=="categorical"){
				PRED <- predict(model.obj, v[!nonPredict,,drop=FALSE])
				
				if(any(is.na(suppressWarnings(as.numeric(Ylev))))){
					v.pred[!nonPredict] <- as.integer(PRED)
				}else{
					v.pred[!nonPredict] <- as.integer(as.character(PRED))
				}
			}

			if(response.type=="continuous"){
				v.pred[!nonPredict] <- predict(model.obj, v[!nonPredict,,drop=FALSE])}
		}
		if(model.type=="CF"){
			if(response.type=="binary"){
				v.pred[!nonPredict] <- signif(CF.list2df(predict(model.obj, newdata=v[!nonPredict,,drop=FALSE],OOB=FALSE,type="prob"))[,2],2)
			}
			if(response.type=="categorical"){
				PRED <- predict(model.obj, newdata=v[!nonPredict,,drop=FALSE], OOB=FALSE, type="response")

				if(any(is.na(suppressWarnings(as.numeric(Ylev))))){
					v.pred[!nonPredict] <- as.integer(PRED)
				}else{
					v.pred[!nonPredict] <- as.integer(as.character(PRED))
				}

			}
			if(response.type=="continuous"){
				v.pred[!nonPredict] <- predict(model.obj, newdata=v[!nonPredict,,drop=FALSE], OOB=FALSE, type="response")
			}
		}
		if(model.type=="QRF"){
			#v.pred[!nonPredict,,drop=FALSE] <- ###doesn't work!
			v.pred[!nonPredict,] <- predict(model.obj.QRF, v[!nonPredict,,drop=FALSE], quantiles=quantiles)
			if(!is.null(model.obj.RF)){
				v.pred.RF[!nonPredict] <- predict(model.obj.RF, newdata=v[!nonPredict,,drop=FALSE], OOB=FALSE, type="response")
				writeValues(out.RF, v.pred.RF, r)
				}
		}
		if(model.type=="SGB"){
			if(response.type=="categorical"){
				 vote <- gbm::predict.gbm(	object=model.obj,
									newdata=v[!nonPredict,,drop=FALSE],
									n.trees=n.trees,
									type="response",
									single.tree=FALSE)[,,1]
				PRED<-factor(colnames(vote)[apply(vote,1,which.max)],levels=Ylev)
				if(any(is.na(suppressWarnings(as.numeric(Ylev))))){
					v.pred[!nonPredict] <- as.integer(PRED)
				}else{
					v.pred[!nonPredict] <- as.integer(as.character(PRED))
				}
			}else{
				v.pred[!nonPredict] <- gbm::predict.gbm(	object=model.obj,
											newdata=v[!nonPredict,,drop=FALSE],
											n.trees=n.trees,
											type="response",
											single.tree=FALSE)
			}
		}

	}

	writeValues(out, v.pred, r)
	
	if(map.sd){
		v.mean  <- rep(NA,length=nrow(v))
		v.stdev <- rep(NA,length=nrow(v))
		v.coefv <- rep(NA,length=nrow(v))

		if(any(!nonPredict)){
			if(model.type=="RF"){
				v.everytree <- predict(model.obj, v[!nonPredict,], predict.all=TRUE)$individual} #only has rows for !nonPredict
			if(model.type=="QRF"){
				v.everytree <- predict(model.obj.RF, v[!nonPredict,], predict.all=TRUE)$individual}
			v.mean[!nonPredict]  <- apply(v.everytree,1,mean)
			v.stdev[!nonPredict] <- apply(v.everytree,1,sd)
			v.coefv[!nonPredict] <- v.stdev[!nonPredict]/v.mean[!nonPredict]
			v.coefv[!nonPredict][v.stdev[!nonPredict]==0]<-0
			v.coefv[!nonPredict][v.mean[!nonPredict]==0]<-0
			rm(v.everytree)
		}

		writeValues(out.mean, v.mean, r)
		writeValues(out.stdev, v.stdev, r)
		writeValues(out.coefv, v.coefv, r)
	}
}

################################################################################################

out <- writeStop(out)

if(model.type%in%c("RF","CF","SGB")){
	out<-setMinMax(out)
}
if(model.type=="QRF"){
	for(i in 1:length(quantiles)){
		out[[i]] <- setMinMax(out[[i]])}
	if(!is.null(model.obj.RF)){
		out.RF<-writeStop(out.RF)
		out.RF<-setMinMax(out.RF)
		writeRaster(out.RF,OUTPUTfn.RF,overwrite=TRUE,datatype=data.type)}
}

writeRaster(out,OUTPUTfn,overwrite=TRUE, datatype=data.type)

if(map.sd){
	out.mean  <- writeStop(out.mean)
	out.stdev <- writeStop(out.stdev)
	out.coefv <- writeStop(out.coefv)
	
	out.mean  <- setMinMax(out.mean)
	writeRaster(out.mean, OUTPUTfn.mean,overwrite=TRUE,datatype=data.type)

	out.stdev <- setMinMax(out.stdev)
	writeRaster(out.stdev,OUTPUTfn.stdev,overwrite=TRUE,datatype=data.type)

	out.coefv <- setMinMax(out.coefv)
	writeRaster(out.coefv,OUTPUTfn.coefv,overwrite=TRUE,datatype=data.type)
}

###clean up tmp files

#print(OUTPUTfn.brick)
#print(TMPfn.map)

file.remove(TMPfn.map)
file.remove(grd2gri(TMPfn.map))

if(!keep.predictor.brick){
	file.remove(OUTPUTfn.brick)
	file.remove(grd2gri(OUTPUTfn.brick))
}
if(model.type=="QRF"){
	if(!is.null(model.obj.RF)){
		file.remove(TMPfn.RF)
		file.remove(grd2gri(TMPfn.RF))
	}
}
if(map.sd){
	file.remove(TMPfn.mean)
	file.remove(TMPfn.stdev)
	file.remove(TMPfn.coefv)
	file.remove(grd2gri(TMPfn.mean))
	file.remove(grd2gri(TMPfn.stdev))
	file.remove(grd2gri(TMPfn.coefv))
}

}

