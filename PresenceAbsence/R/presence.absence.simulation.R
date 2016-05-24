
presence.absence.simulation<-function(	n,
							prevalence,
							N.models=1,
							shape1.absent=1,shape2.absent=1, 
							shape1.present=1,shape2.present=1){
###Simulates presence-absence data. 
###First, Observed values are generated with a binomial function,
###then two beta distributions are used to generate predicted values.
###
###	Shape arguments give the parameters for the beta distributions.  
###	If N.models=1 the shape arguments should be single numbers. 
###	If N.models>1, the shape arguments should be vectors of length N.modles.
###
###	The mean of the beta distribution equals shape1/(shape1+shape2)
###	To get reasonable predictions (e.g. better than random), the mean for the plots where the
###	observed value is present should be higher than that of the plots where the species is absent.
###
###	n			number of plots
###	prevalence		probability of speciece presence for observed values
###	N.models		number of models to simulate predictions for
###	shape1.absent	first parameter for beta distribution for plots where observed value is absent
###	shape2.absent	second parameter for beta distribution for plots where observed value is absent
###	shape1.present	first parameter for beta distribution for plots where observed value is present
###	shape2.present	second parameter for beta distribution for plots where observed value is present

###Check that N.models is valid

if(N.models!= round(N.models) || N.models<1){
	stop("N.models must be a positive integer")}

###Check that 'prevalence' is valid

if(prevalence<0 || prevalence>1){
	stop("'prevalence' must be a probability between zero and one")}

###Check if shape parameters are valid

if(min(shape1.absent,shape2.absent,shape1.present,shape1.present)<=0){
	stop("shape parameters must be positive")}

###Check length of parameters equals number of models

if(length(shape1.absent)!= N.models){
	if(length(shape1.absent)!=1){
		stop("shape parameters must of length 1 or length equal to number of models to generate")
	}else{shape1.absent<-rep(shape1.absent,N.models)}
}
		
if(length(shape2.absent)!= N.models){
	if(length(shape2.absent)!=1){
		stop("shape parameters must of length 1 or length equal to number of models to generate")
	}else{shape2.absent<-rep(shape2.absent,N.models)}
}

if(length(shape1.present)!= N.models){
	if(length(shape1.present)!=1){
		stop("shape parameters must of length 1 or length equal to number of models to generate")
	}else{shape1.present<-rep(shape1.present,N.models)}
}

if(length(shape2.present)!= N.models){
	if(length(shape2.present)!=1){
		stop("shape parameters must of length 1 or length equal to number of models to generate")
	}else{shape2.present<-rep(shape2.present,N.models)}
}
###Set up DATA

DATA<-as.data.frame(matrix(0,n,2+N.models))
names(DATA)<-c("plotID","Observed",paste("Predicted",1:N.models,sep=""))
DATA$plotID<-1:n

###Simulate Observed values

DATA$Observed<-rbinom(n=n,size=1,prob=prevalence)

###Simulate Predicted values

for(i in 1:N.models){
	DATA[DATA$Observed==1,i+2]<-rbeta(	n=length(DATA[DATA$Observed==1,i+2]),
							shape1=shape1.present[i],
							shape2=shape2.present[i])
	DATA[DATA$Observed==0,i+2]<-rbeta(	n=length(DATA[DATA$Observed==0,i+2]),
							shape1=shape1.absent[i],
							shape2=shape2.absent[i])
	}
return(DATA)
}	

