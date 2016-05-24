#Data must have rows of length Hz*RmL*0.5 on each end that will be trimmed in the final output.
DeadReckoning <- function(rawdata, betas, decinc, Hz=16, RmL=2, DepthHz=1, SpdCalc=1, MaxSpd=NULL)#RmL=runningmean length
#SpdCalc=1 given at same Hz at Mag and Acc data, =2 Given at lower Hz than mag and acc data, =3 calculate via surge integration, =4 Calculate a constant speed, =5 Calculate via depth ***not implemented yet***
{
	#timer=proc.time()[3] # used to determine slow parts of the program and show some output while it is running to indicate where a failure occurs.
	
	#eliminate data where accelerometer data are collected at a speed faster than magnetometer data
	rawdata<-rawdata[!is.na(rawdata$MagSurge),]
	
	Declination=-1*decinc[1]/360*2*pi # formats to Radians
	Inclination=1*decinc[2]/360*2*pi # formats to Radians
	firstdepth<-0 #Initialize the first depth
	
	dims=dim(rawdata)# dimensions need to create fillable matricies
	standardize1=matrix(0,dims[1],6)
	#Standardize the mag and acc channels 
	standardize1[,1]=betas["B0 Intercept","MagSurge"]+betas["B1 Slope","MagSurge"]*rawdata$MagSurge #betas[1,1]+betas[2,1]
	standardize1[,2]=betas["B0 Intercept","MagHeave"]+betas["B1 Slope","MagHeave"]*rawdata$MagHeave #betas[1,2]+betas[2,2]
	standardize1[,3]=betas["B0 Intercept","MagSway"] +betas["B1 Slope","MagSway"] *rawdata$MagSway  #betas[1,3]+betas[2,3]
	standardize1[,4]=betas["B0 Intercept","AccSurge"]+betas["B1 Slope","AccSurge"]*rawdata$AccSurge #betas[1,4]betas[2,4]
	standardize1[,5]=betas["B0 Intercept","AccHeave"]+betas["B1 Slope","AccHeave"]*rawdata$AccHeave #betas[1,5]+betas[2,5]
	standardize1[,6]=betas["B0 Intercept","AccSway"] +betas["B1 Slope","AccSway"] *rawdata$AccSway  #betas[1,6]+betas[2,6]
	#print(c("1",proc.time()[3]-timer)) # Testing and failure information
	 
	# Calculate the running mean over runningmean length (two seconds of data) for the acc data
	runningmean=matrix(0,dims[1],6)
	runningmean[(Hz*(0.5*RmL)):(dims[1]-Hz*(0.5*RmL)),4]=filter(standardize1[1:dims[1],4],rep(1,(RmL*Hz*(0.5*RmL)))/(RmL*Hz*(0.5*RmL)))[(Hz*(0.5*RmL)):(dims[1]-Hz*(0.5*RmL))]
	runningmean[(Hz*(0.5*RmL)):(dims[1]-Hz*(0.5*RmL)),5]=filter(standardize1[1:dims[1],5],rep(1,(RmL*Hz*(0.5*RmL)))/(RmL*Hz*(0.5*RmL)))[(Hz*(0.5*RmL)):(dims[1]-Hz*(0.5*RmL))]
	runningmean[(Hz*(0.5*RmL)):(dims[1]-Hz*(0.5*RmL)),6]=filter(standardize1[1:dims[1],6],rep(1,(RmL*Hz*(0.5*RmL)))/(RmL*Hz*(0.5*RmL)))[(Hz*(0.5*RmL)):(dims[1]-Hz*(0.5*RmL))]
	runningmean[(Hz*(0.5*RmL)):(dims[1]-Hz*(0.5*RmL)),1:3]=standardize1[(Hz*(0.5*RmL)):(dims[1]-Hz*(0.5*RmL)),1:3] # Don't need column means for mag data but don't want to change the rest of the program
	#print(c("2",proc.time()[3]-timer)) # Testing and failure information
		
	dynacc=matrix(0,dims[1],3)
	dynacc=standardize1[,4:6]-runningmean[,4:6] # Calculate the dynamic acceleration
	runningmean[,4:6]=ifelse(runningmean[,4:6]>1,1,ifelse(runningmean[,4:6]<(-1),-1,runningmean[,4:6]))
	rm(standardize1)
	
	####SPEED CALCULATIONS####
	
	# Calculate the x dimension of overall dynamic body acceleration to get a proxy for speed
	xaccdba=matrix(0,dims[1],1)
	xaccdba[Hz:(dims[1]-Hz)]=cumsum(dynacc[Hz:(dims[1]-Hz),1])
	#ODBA=matrix(0,dims[1],1)
	#ODBA=rowSums(abs(dynacc))
	rm(dynacc)
	
	# Standardize X odba acceleration to between 0 and 3.5 meters per sec
	minxacc=min(xaccdba)
	maxxacc=max(xaccdba)
	
	#if speed is given
	if(SpdCalc==1)
	{
		calcspd<-rawdata$Speed
	}
	#if speed is given at a Hz lower than Magnetometers
	if(SpdCalc==2)
	{
		calcspd<-approx(seq(1:dim(rawdata)[1]),rawdata$Speed,xout=seq(1:dim(rawdata)[1]),rule=2)$y
	}
	#For a linear function of speed to odba use....
	if(SpdCalc==3)
	{
		speedmat<-matrix(data= c(minxacc,maxxacc,0,MaxSpd),nrow=2,ncol=2,byrow=FALSE)
		spdmdl<-lm(speedmat[,2]~speedmat[,1])
		calcspd<-matrix(0,dims[1],1)
		calcspd[(Hz*(0.5*RmL)):(dims[1]-Hz*(0.5*RmL)),1]=(xaccdba[(Hz*(0.5*RmL)):(dims[1]-Hz*(0.5*RmL))]*spdmdl$coef[2]+spdmdl$coef[1])/Hz
		calcspd[(Hz*(0.5*RmL)):(dims[1]-Hz*(0.5*RmL)),1]=calcspd[(Hz*(0.5*RmL)):(dims[1]-Hz*(0.5*RmL)),1]+match(calcspd[(Hz*(0.5*RmL)):(dims[1]-Hz*(0.5*RmL)),1],0,nomatch=0)*0.0000001 #This turns speeds of 0 into 0.0000001
	}
	#For a power function of speed to odba use...
	#maxspeed=3.5
	#spd<- function(x)
	#{
	#	((minxacc+maxspeed^x)-maxxacc)^2
	#}
	#spdmdl<-optimize(spd, c(1,4), tol=0.0001)
	#calcspd=matrix(0,dims[1],1)
	#calcspd[16:(dims[1]-16),1]=((xaccdba[16:(dims[1]-16)]-minxacc)^(1/spdmdl$minimum))/16 #!!!!BUT CHECK THIS LINE
	#calcspd[16:(dims[1]-16),1]=calcspd[16:(dims[1]-16),1]+match(calcspd[16:(dims[1]-16),1],0,nomatch=0)*0.0000001 #This turns speeds of 0 into 0.0000001
	#Assume constant speed
	if(SpdCalc==4)
	{
		calcspd<-matrix(MaxSpd,dims[1],1)
		
		
	}
	#Resting Behavior given a speed of 0.000002.
	#behavior=matrix(0,dims[1],1)
	#DateTime=paste(as.character(rawdata[,1]),as.character(rawdata[,2])) #Combine date and time
	#matchnum=match(DateTime[1],RDBeh[,1])+3
	#behavior[16:(dims[1]-16),1]<-RDBeh[(matchnum+16):(matchnum+(dims[1]-16)),2]
	#calcspd[which(behavior==1)]<-0.000002
	#rm(behavior)
	#rm(DateTime)
	
	#print(c("4",proc.time()[3]-timer)) # Testing and failure information
	
	rm(xaccdba)
	
	#Fill in the depth between seconds as a linear interpolation
	#Tiph approx
	calcdpth<-approx(seq(1:dim(rawdata)[1]),rawdata$Depth,xout=seq(1:dim(rawdata)[1]),rule=2)$y
	
	#print(c("5",proc.time()[3]-timer)) # Testing and failure information
	
	# DeadReckoning
	hvec=c(cos(pi/2-Declination)*cos(Inclination),sin(pi/2-Declination)*cos(Inclination),-sin(Inclination))

	gvec=c(0,0,-1)
	hxgvec=c(hvec[2]*gvec[3]-hvec[3]*gvec[2],hvec[3]*gvec[1]-hvec[1]*gvec[3],hvec[1]*gvec[2]-hvec[2]*gvec[1])
	hg=sum(hvec*gvec)
	#print(c("5.1",proc.time()[3]-timer)) # Testing and failure information
	
	runningmean[,4:6]=ifelse(runningmean[,4:6]>1,1,ifelse(runningmean[,4:6]<(-1),-1,runningmean[,4:6]))
	
	xmat=matrix(0,dims[1],3)
	xmat[,1]=(runningmean[,1]-runningmean[,4]*hg)/(1-hg^2)
	xmat[,2]=(runningmean[,4]-runningmean[,1]*hg)/(1-hg^2)
	ymat=matrix(0,dims[1],3)
	ymat[,1]=(runningmean[,2]-runningmean[,5]*hg)/(1-hg^2)
	ymat[,2]=(runningmean[,5]-runningmean[,2]*hg)/(1-hg^2)
	zmat=matrix(0,dims[1],3)
	zmat[,1]=(runningmean[,3]-runningmean[,6]*hg)/(1-hg^2)
	zmat[,2]=(runningmean[,6]-runningmean[,3]*hg)/(1-hg^2)
	rm(runningmean)
	
	#print(c("5.5",proc.time()[3]-timer)) # Testing and failure information
	xmat[,3]=ymat[,1]*zmat[,2]-ymat[,2]*zmat[,1] #old XandY
	#xmat[,3]=zmat[,1]*ymat[,2]-zmat[,2]*ymat[,1] #new XandY
	rm(ymat)
	rm(zmat)
	#print(c("5.6",proc.time()[3]-timer)) # Testing and failure information
	xmat2=matrix(0,dims[1],3)
	xmat2[,1]=xmat[,1]*hvec[1]+xmat[,2]*gvec[1]+xmat[,3]*hxgvec[1]
	xmat2[,2]=xmat[,1]*hvec[2]+xmat[,2]*gvec[2]+xmat[,3]*hxgvec[2]
	xmat2[,3]=xmat[,1]*hvec[3]+xmat[,2]*gvec[3]+xmat[,3]*hxgvec[3]
	rm(xmat)
	#rm(gvec)
	rm(hvec)
	rm(hxgvec)
	rm(hg)
	
	headingrad=(atan2(xmat2[,1],xmat2[,2])+ifelse(atan2(xmat2[,1],xmat2[,2])<0,2*pi,0)) # Heading in radians
	
	rm(xmat2)
	x1=sin(headingrad)*calcspd
	y1=cos(headingrad)*calcspd
	rm(headingrad)
	#Make output matrix for xyz dimension data
	#print(c("8",proc.time()[3]-timer)) # Testing and failure information
	xyz=matrix(0,dims[1],5,dimnames=list(1:dims[1],c("DateTime","Xdim","Ydim","Depth","Speed")))
	xyz[,4]=calcdpth
	rm(calcdpth)
	#xyz[,6]=ODBA
	#rm(ODBA)
	xyz[,5]=calcspd
	rm(calcspd)
		
	xyz[,3]=cumsum(y1)
	xyz[,2]=cumsum(x1)
	rm(x1)
	rm(y1)
	
	#print(c("9",proc.time()[3]-timer)) # Testing and failure information
	xyz=as.data.frame(xyz,colclasses=c(DateTime="character",Xdim="double",Ydim="double",Depth="double",ODBA="double"))
	#print(c("9.1",proc.time()[3]-timer)) # Testing and failure information
	if('DateTime' %in% colnames(rawdata))
	{
		xyz[,1]<-rawdata$DateTime
	}else{
		xyz[,1]<-paste(as.character(rawdata$Date),as.character(rawdata$Time)) #Combine date and time
	}
	#print(c("9.2",proc.time()[3]-timer)) # Testing and failure information
	xyz=xyz[(Hz-1):(dims[1]-Hz),]
	
	#print(c("10",proc.time()[3]-timer)) # Testing and failure information	
	
	return(xyz=xyz)
	
}

