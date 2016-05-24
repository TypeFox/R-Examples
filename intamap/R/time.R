#estimation of time model 
estimateTimeModel<-function(FUN,class,formulaString,debug.level = 0,...){
	dots = list(...)
	times<-timeSampling(FUN=FUN,class=class,formulaString=formulaString,...)
	if (debug.level >= 1) {
    write.table(times$eTimes,file = paste(FUN,class,".txt",sep=""),append=TRUE)
  	write.table(times$pTimes,file = paste(FUN,class,".txt",sep=""),append=TRUE)
	}
  eTimeModel=loess(data=data.frame(times$eTimes),Time~nObs,surface="direct")
  if (class == "copula") {
    pTimeModel=loess(data=data.frame(times$pTimes),Time~nObs+nPred+nType,surface="direct")
	} else pTimeModel=loess(data=data.frame(times$pTimes),Time~nObs+nPred,surface="direct")
  return(list(eTimeModel = eTimeModel, pTimeModel = pTimeModel, times = times))
}

#function that gives the time estimation
predictTime=function(nObs,nPred,class,formulaString,calibration=FALSE,outputWhat,FUN = "spatialPredict",...){
	if (exists("timeModels", envir = .GlobalEnv)){ 
		cat("[Time models loaded...]\n")
  	timeModels = get("timeModels", envir = .GlobalEnv)
	} else { 
#			rpath=system.file("models",package="intamap")
#			timeModels=try(load(paste(rpath,"/timeModel.Rdata",sep="")))
#			if(inherits(timeModels,"try-error")) 	timeModels=NULL
#			else	
#  	timeModels = list()
    data(timeModels, envir = environment())
    warning(paste("\n using standard model for estimating time. For better \n",
         "platform spesific predictions, please run \n",
         "timeModels <- generateTimeModels()\n ",
         "and save the workspace"))
	}
	
	#load dots in case of spatialPredict function
	dots=list(...)
	
	if(FUN=="spatialPredict"){ 
		if(class %in% c("linearVariogram","automap")){
				if(formulaString==as.formula(value~1)) FUN=paste(FUN,class,"OK",sep="_")
				else	FUN=paste(FUN,class,"UK",sep="_")
				}else{
					FUN=paste(FUN,class,sep="_")
			}
	}

 	eModels=timeModels[[as.character(FUN)]]$eTimeModel
	if (!is.null(eModels)) {
	  tmm = as.data.frame(timeModels[[as.character(FUN)]]$times$pTimes)
    nObsMax = max(tmm$nObs[tmm$nPred > 0])
    nPredMax = max(timeModels[[as.character(FUN)]]$times$pTimes[,"nPred"])
    eTime = c(1:2)
    pTime = c(1:2)
    for (ite in 1:2) {
      if (ite == 2 && nObs > nObsMax) nObs = nObsMax
      if (ite == 2 && nPred > nPredMax) nObs = nPredMax
      eTime[ite]<-predict(eModels,data.frame(nObs=nObs),se=TRUE)[[1]]
      if (class == "copula") {
      	nType = length(names(outputWhat)[!(names(outputWhat) %in% c("mean","variance"))])
        pModels=timeModels[[as.character(FUN)]]$pTimeModel
	      pTime[ite]<-predict(pModels,data.frame(nObs=nObs,nPred=nPred,nType = nType),se=TRUE)[[1]] 	  
     	} else {
        pModels = timeModels[[as.character(FUN)]]$pTimeModel
      	pTime[ite]<-predict(pModels,data.frame(nObs=nObs,nPred=nPred),se=TRUE)[[1]]
      }
    }
    
    predTime = max(eTime) + max(pTime)
    if (predTime < 2) predTime = 2
#  ind=which(result$fit<2)
#	result$fit[ind]=2
  } else {
    warning(paste("Could not find time model for method",FUN))
    predTime = NA
  }
  return(predTime)
#	if(calibration) return ((result$fit)/(timeModels$cl)*timeModels$clnew)
#		else	return(result$fit)
}




#function that does machine calibration
timeCalibration = function(save=FALSE,original=TRUE){
	if(exists("timeModels")){ 
		get("timeModels")
		cat("[Time models loaded...]\n")
	}else{ 
			rpath=system.file("data",package="intamap")
#			timeModels=try(load(paste(rpath,"/timeModel.Rdata",sep="")))
#			if(inherits(timeModels,"try-error")) 	
			timeModels=NULL
			}
	
	load(file=rpath)
	nIter=100
	n=50

  grd=expand.grid(x=seq(1,n),y=seq(1,n))
  z=rnorm(1:n^2)
  dataSim=data.frame(x=grd$x,y=grd$y,z)
  mes=dataSim[sample(1:n^2,n),]
  coordinates(grd)=~x+y
  coordinates(dataSim)=~x+y
  gridded(grd)=T

  t1=as.numeric(proc.time()[3])
  for (i in 1:nIter) k=krige(z~1,dataSim,grd)
  cl=as.numeric((proc.time()[3])-t1)/nIter

	if(save){
	
		if(original) timeModels$cl=cl
			else timeModels$clnew=cl
	rm(list=ls(all.names=T)[-which((ls(all.names=T))=="timeModels")])
	save.image()
	}
	return(cl)
}  

############################################################
#Generate the ellipse
############################################################
ellipse2 = function (x, scale = c(1, 1), centre = c(0, 0), level = 0.95,
   t = sqrt(qchisq(level, 2)), which = c(1, 2), npoints = 100,
    ...)
{
    names <- c("x", "y")
    if (is.matrix(x)) {
        xind <- which[1]
        yind <- which[2]
        r <- x[xind, yind]
        if (missing(scale)) {
            scale <- sqrt(c(x[xind, xind], x[yind, yind]))
            if (scale[1] > 0)
                r <- r/scale[1]
            if (scale[2] > 0)
                r <- r/scale[2]
        }
        if (!is.null(dimnames(x)[[1]]))
            names <- dimnames(x)[[1]][c(xind, yind)]
    }
    else r <- x
    r <- min(max(r, -1), 1)

    dd <- acos(0)
   # print(dd)
    a <- seq(0, 2 * pi, len = npoints)
    unrot = matrix(c(t * scale[1] * cos(a + dd/2) + centre[1], t * scale[2] *
        cos(a - dd/2) + centre[2]), npoints, 2, dimnames = list(NULL,
        names))

   rotm = matrix(c(cos(r),-sin(r),sin(r),cos(r)),ncol=2)
   #print(rotm)
   #print(unrot)
   rotmat = unrot %*% rotm
   #print(rotmat)
   rotmat
}

###################################################
#function generates random intamap object
###################################################
generateRandomObject<-function( stot = 50, nPred=1000,FUN,class,formulaString,...){
	# stot = 300
	nclus = 5
 	cnum = as.integer(stot/(nclus+1))
 	boun = data.frame(x=c(0,1000,1000,0,0),y=c(0,0,1000,1000,0))
 	coordinates(boun) = ~x+y
 	boun = Polygon(boun)
 	clus=spsample(boun,nclus,"random")
 	predictionLocations = spsample(boun,cnum,"random")
 for (i in 1:5) {
   ratio = 10
   while (ratio>3) {
     clsize1 = runif(1,0,150)
     clsize2 = runif(1,0,150)
     ratio = clsize2/clsize1
     if (ratio<1) ratio = 1/ratio
     ang = runif(1,-1,1)
	  ellip = as.data.frame(ellipse2(ang,scale = c(x=clsize1,y=clsize2),centre=coordinates(clus[i,]),npoints=20))
     names(ellip) = c("x","y")
     ellip = ellip[ellip$x > 0 & ellip$y > 0 & ellip$x < 1000 & ellip$y<1000,]
  	   if (any(is.na(ellip)) | dim(ellip)[1]<5) ratio = 100
  	 }
  	 ellip = rbind(ellip,ellip[1,])
  	 #plot(ellip)
  	 coordinates(ellip) = ~x+y
  	 ellip = Polygon(ellip)
  	 elsamp = spsample(ellip,cnum,"random")
  	 #plot(elsamp,asp=1)
  	 predictionLocations = rbind(predictionLocations,elsamp)
	}
	  if (dim(zerodist(predictionLocations,zero = 0))[1] > 0)
     predictionLocations = remove.duplicates(predictionLocations,zero = 1.0)		


		p1 = SpatialPointsDataFrame(predictionLocations[sample(dim(coordinates(predictionLocations))[1],1)],data.frame(value = 0))
 		observations = krige(value~1,p1,predictionLocations,vgm(1,"Ste",500,0.3,kappa = 0.3),nsim=1,nmax=50)
   	names(observations) = "value"
 		#generate intamap Object with ellipse data
		if (class == "transGaussian") observations$value = exp(observations$value)
  	obj= createIntamapObject(observations = observations,
           predictionLocations = spsample(boun,nPred,"random"),
           class=class, formulaString=formulaString, ...)
	



		#generate intamap Object with ellipse data
		#	len=dim(coordinates(predictionLocations))[1] 	
		#	obj= createIntamapObject(observations = SpatialPointsDataFrame(coordinates(predictionLocations),data.frame(value=runif(len))),
  	#																			predictionLocations = spsample(boun,nPred,"random"),formula=formula,class=class,...)



#		if(FUN=="spatialPredict"){
	
#		if(class=="copula"){
#				print("Copula")
#				obj=estimateParameters(obj)
				#the worst case scenario for copulas !
#			obj$params$copulaParams$margin<- list(name="norm",params=3)
#			obj$params$copulaParams$correlation<-list(method="Ste", params=c(0.01,30,3))
#			obj$params$copulaParams$trend<-list(F=as.matrix(1,200),params=20)
#			obj=estimateAnisotropy(obj)
	
#				obj$params$copulaParams$anisotropy<-list(params=c(obj$anisPar$direction*pi/180,obj$anisPar$ratio))
#				obj$params$copulaParams$copula<-list(method="norm") #another option would be "chisq"
#		} 
#	}


return(obj)
}

#internal function for time sampling
timeSampling<-function(FUN,class,formulaString, nSam = 1,...){

  if (class == "copula") nSam = max(nSam %/%2,1)	
  obsInterv=setTimeObsInterv(FUN,class, ...)
	predInterv=setTimePredInterv(FUN,class, ...)

  compNum=length(obsInterv)*length(predInterv) *nSam
	
	timing<-matrix(0,compNum,3)
	print("Timing Started")	
              
	outputWhat0 = list(mean = TRUE, variance = TRUE, quantile = 0.5, excprob = 3, excprob = 0, excprob = -1)
  if (class == "idw") outputWhat0 = list(mean = TRUE)
  if (class == "copula") nOutput = length(outputWhat0) else nOutput = 2
  eTimes = matrix(0,length(obsInterv)*nSam,2)
  if (class == "copula") {
    compNum = compNum*(nOutput-1)
    nCol = 4
  } else nCol = 3
  pTimes = matrix(0,compNum,nCol)
	ie = 0
	ip = 0
	for (ii in 1:nSam) {
    for (i in 1:(length(obsInterv)-ifelse(class == "copula",2,0))){
  		ie = ie + 1
      obj<-generateRandomObject(obsInterv[i],predInterv[1],FUN=FUN,
          formulaString=formulaString,class=class,...)
      eTime=system.time(obj <- estimateParameters(obj))
      eTimes[ie,1] = obsInterv[i]
      eTimes[ie,2] = eTime[1]
      for (j in 1:length(predInterv)){
        obj$predictionLocations = spsample(obj$observations,predInterv[j],"regular")
#			obj<-generateRandomObject(obsInterv[i],predInterv[j],FUN=FUN,formula=formula,class=class,...)
#			if (!is.null(obj$params$nmax)) nmax=obj$params$nmax else
# 	    nmax=dim(obj$observations)[1]-1
        for (k in 2:nOutput) {
          ip = ip + 1
          if (class == "copula") outputWhat = outputWhat0[1:k] else outputWhat = outputWhat0
          obj$outputWhat = outputWhat
          pTime = try(system.time(obj <- spatialPredict(obj)))
          if (!is(pTime,"try-error")) {
	      		print(pTime)
            pTimes[ip,1] = obsInterv[i]
            pTimes[ip,2] = predInterv[j]
            pTimes[ip,nCol] = pTime[1]
            if (class == "copula") pTimes[ip,3] = max(0,k-2)
            print(paste("Found time for FUN = ",FUN,"class",class,
                    "iSam",ii,"of",nSam,"obs",i,"of ",length(obsInterv),
                    "pred", j,"of ",length(predInterv), "nOoutput:",k))
          }
        }
  		}	
  	}
  }
  pTimes = pTimes[pTimes[,1] > 0,]
  eTimes0 = matrix(c(0,0),ncol = 2)
  if (class == "copula") {
    for (i in 1:max(length(predInterv),length(obsInterv))) {
      for (j in 0:2) {
        if (i == 1 & j == 0) {
          pTimes0 = matrix(c(obsInterv[i],0,j,0),nrow=1) 
        } else if (i <= length(obsInterv)) {
          pTimes0 = rbind(pTimes0,c(obsInterv[i],0,j,0))
        } 
        if (i <= length(predInterv)) pTimes0 = rbind(pTimes0,c(0,predInterv[i],j,0))
      }
    }          
  	columnNames=c("nObs","nPred","nType","Time")
  } else {
  	pTimes0 = matrix(c(obsInterv,
                     rep(0,length(predInterv)),
                     rep(0,length(obsInterv)),
                     predInterv,
                     rep(0,length(obsInterv)+length(predInterv))),ncol = 3)
  	columnNames=c("nObs","nPred","Time")
  }
  eTimes = rbind(eTimes0,eTimes)
  pTimes = rbind(pTimes0,pTimes)

	dimnames(pTimes)[[2]]=columnNames	
	dimnames(eTimes)[[2]]=columnNames[c(1,nCol)]		
	return(list(eTimes = eTimes, pTimes = pTimes))
}

#internal function setting the observations intervals
setTimeObsInterv=function(FUN,class, test=FALSE, ...){
	obs=switch(FUN,
				"estimateAnisotropy"=10*c(10,50,100,200,500),
				"spatialPredict"=switch(class, "copula"=c(20,50,100,200,400, 800),
																	c(20,100,300,1000,2000,5000,10000)),
				"doSegmentation"=seq(200,5000,by=500),
				c(25,50,100,200,500)
			)
  if (test) obs = obs[1:4]/8
	return(obs)
}
#internal function setting the prediction locations intervals
setTimePredInterv=function(FUN,class, test=FALSE, ...){
	pred=switch(FUN,
				"estimateAnisotropy"=c(1,500,1000,2000,5000,10000),
				"spatialPredict"=c(100,500,1000,2000,5000,10000, 20000, 50000),
				"doSegmentation"=c(1,500,1000,2000,5000,10000),
				c(10,500,1000,2000,5000,10000)
			)

	if (test) pred = pred[1:4]/8
  return(pred)
}



#function that iteratively estimates the time models 
generateTimeModels=function(genClasses = NULL, noGenClasses = NULL,nSam = 1, test = FALSE, debug.level = 0){
#	require(ellipse)
#	namespace=c("estimateAnisotropy","doSegmentation")
#	formulas=c("value~1","value~r1+r2")
	formulas=c("value~1")
	classes=c("linearVariogram","automap")
	classes2=c("transGaussian","copula","psgp","idw")
	timeModels=NULL
  if (!is.null(genClasses)) {
    classes = classes[classes %in% genClasses]
    classes2 = classes2[classes2 %in% genClasses]
  }
  if (!is.null(noGenClasses)) {
    classes = classes[!(classes %in% noGenClasses)]
    classes2 = classes2[!(classes2 %in% noGenClasses)]
  }
	#calculate the rest time models 
#	for(l in 1:length(namespace)){
#		print(l)
#		timeModels[[namespace[l]]]=estimateTimeModel(FUN=namespace[l],
#																								 formula=value~1,
#																								 class="automap")
#	}

	#calculate Models for spatialPredict function
	for(j in 1:length(classes)){
		for (i in 1:length(formulas)){

			print(paste(i,j,sep="-"))

			if(i==1)	model=paste("spatialPredict",classes[j],"OK",sep="_")	
			if(i==2)		model=paste("spatialPredict",classes[j],"UK",sep="_")	

			if (j > 0 && !is.na(classes[j])) 
          timeModels[[model]] = estimateTimeModel(FUN="spatialPredict",
               formulaString=as.formula(formulas[i]), nSam = nSam, 
               class=classes[j], test = test, debug.level = debug.level )
		}
	}


	for(k in 1:length(classes2)){
			model=paste("spatialPredict",classes2[k],sep="_")	
			if (k >0 && !is.na(classes2[k])) 
              timeModels[[model]] = estimateTimeModel(FUN="spatialPredict",
              formulaString="value~1", nSam = nSam, class=classes2[k], 
              test = test, debug.level = debug.level )

	}

	return(timeModels)
}	
