Standardize <-
function(MagOrSR,MagOrHV,MagOrSW,AccOrSR,AccOrHV,AccOrSW,magSRmin,magSRmax,magHVmin,magHVmax,magSWmin,magSWmax,accSRmin,accSRmax,accHVmin,accHVmax,accSWmin,accSWmax)
{
	#Standardization/Calibrating/normalizaiton/relativizing acceleration and magnetic data
	#The right-hand-rule indicates the polarity of the magnetometers and accelerometers required for the pseudotrack reconstruction algorithm, so that when the front, top or left side 
	#of the tag is facing the earth, the accelerometers are at their minimal (-1) reading.  A similar rule applies to the magnetic field except that the maximal (+1) reading of the tag 
	#will be facing north and at the angle of inclination of the magnetic field (the angle at which the magnetic field enters the earth) at that place on the planet.  If the tag sensors
	#conform to this rule, the *Or* parameters should all be 1.  If any sensors do not conform to this, they appropriate *Or* parameter should be -1.  Instead of labeling the tag in 
	#X,Y, and Z dimensions, the directions are labeled as Surge, Heave and Sway axes where the Surge indicates the front(anterior) and back(posterior) axis, the Heave is the top(forsal)
	#and bottom (ventral) axis and the Sway is the right and left (lateral) axis.  

	Betas=matrix(0,2,6)
	#Surge Magnetometer standardization
	stdMSRmat=matrix(data= c(magSRmin,magSRmax,MagOrSR*-1,MagOrSR*1),nrow=2,ncol=2,byrow=FALSE)
	stdMSRmdl<-lm(stdMSRmat[,2]~stdMSRmat[,1])
	Betas[1,1]<-stdMSRmdl$coef[1]
	Betas[2,1]<-stdMSRmdl$coef[2]
	#Heave Magnetometer standardization
	stdMHVmat=matrix(data= c(magHVmin,magHVmax,MagOrHV*-1,MagOrHV*1),nrow=2,ncol=2,byrow=FALSE)
	stdMHVmdl<-lm(stdMHVmat[,2]~stdMHVmat[,1])
	Betas[1,2]<-stdMHVmdl$coef[1]
	Betas[2,2]<-stdMHVmdl$coef[2]	
	#Sway Magnetometer standardization
	stdMSWmat=matrix(data= c(magSWmin,magSWmax,MagOrSW*-1,MagOrSW*1),nrow=2,ncol=2,byrow=FALSE)
	stdMSWmdl<-lm(stdMSWmat[,2]~stdMSWmat[,1])
	Betas[1,3]<-stdMSWmdl$coef[1]
	Betas[2,3]<-stdMSWmdl$coef[2]	
	#Surge Accelerometer standardization
	stdASRmat=matrix(data= c(accSRmin,accSRmax,AccOrSR*-1,AccOrSR*1),nrow=2,ncol=2,byrow=FALSE)
	stdASRmdl<-lm(stdASRmat[,2]~stdASRmat[,1])
	Betas[1,4]<-stdASRmdl$coef[1]
	Betas[2,4]<-stdASRmdl$coef[2]	
	#Heave Accelerometer standardization
	stdAHVmat=matrix(data= c(accHVmin,accHVmax,AccOrHV*-1,AccOrHV*1),nrow=2,ncol=2,byrow=FALSE)
	stdAHVmdl<-lm(stdAHVmat[,2]~stdAHVmat[,1])
	Betas[1,5]<-stdAHVmdl$coef[1]
	Betas[2,5]<-stdAHVmdl$coef[2]	
	#Sway Accelerometer standardization
	stdASWmat=matrix(data= c(accSWmin,accSWmax,AccOrSW*-1,AccOrSW*1),nrow=2,ncol=2,byrow=FALSE)
	stdASWmdl<-lm(stdASWmat[,2]~stdASWmat[,1])
	Betas[1,6]<-stdASWmdl$coef[1]
	Betas[2,6]<-stdASWmdl$coef[2]
	colnames(Betas)<-c("MagSurge","MagHeave","MagSway","AccSurge","AccHeave","AccSway")
	rownames(Betas)<-c("B0 Intercept","B1 Slope")

	return=(Betas=Betas)
}
