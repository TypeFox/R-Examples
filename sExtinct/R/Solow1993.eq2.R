Solow1993.eq2 <-
function(sightingdata, alpha, test.year, data.out){

##some warnings

##data frame warnings
   if(missing(sightingdata)) stop('no data present', call.=FALSE)
   
   if(!is.numeric(sightingdata[,1])) 
   		stop('sighting data must be a data.frame of real numbers', call.=FALSE)
   		
   if(!is.numeric(sightingdata[,2]))
      	stop('sighting data must be a dataframe of real numbers', call.=FALSE)
      	
   if(ncol(sightingdata)<2) 
   		stop('some data missing', call.=FALSE)
   
   if(max(sightingdata[,2])>max(sightingdata[,1]))
      warning("Check column values are correct, values for number of sighting events appear to be larger than the years at which sightings occur", call.=FALSE)

   if(ncol(sightingdata)>2)
   		warning("data passed to the function is larger than 2 columns, extinction is being estimated from the values in the first two columns only", call.=FALSE)
   		
   if(min(sightingdata[,2])<0)
	   stop('sightingdata cannot contain negative values', call.=FALSE)   	
	   
	if(min(sightingdata[,1])<0)
	   stop('sightingdata cannot contain negative values', call.=FALSE)  
	   
   if(min(diff(sightingdata[,1]))==0) 
		stop('multiple identical year values entered', call.=FALSE)
		
   if(min(diff(sightingdata[,1]))<0)
		warning("sightingdata has been ordered by the time column", call.=FALSE)
		
   if(min(diff(sightingdata[,1]))<0)
   		sightingdata<-sightingdata[order(sightingdata[,1]),]
   		
##missing value warnings
   if(missing(alpha)) 
   		stop('no alpha value specified', call.=FALSE)
   
   if(alpha>1) 
   		stop('alpha value must be less than 1', call.=FALSE)
   
   if(alpha<0) 
   		stop('alpha value must be a positive number', call.=FALSE)
   
   if(!is.numeric(test.year)) 
   		stop('no "test.year" specified', call.=FALSE)
   
   if(missing(data.out)) 
		stop('data.out must be specified as either TRUE or FALSE', call.=FALSE)
   
   if(!is.logical(data.out)) 
   		stop('data.out must be logical', call.=FALSE)


	
	names(sightingdata)<-c("yrs", "sights")
	sightingdata<-subset(sightingdata, sightingdata$sights>0)
	yrs<-seq(sightingdata[1,1], test.year) 
	sights<-rep(0, times=length(yrs))
	res<-data.frame(yrs, sights)
	for(i in 1:length(sightingdata$yrs)){res$sights[length(seq(sightingdata$yrs[1], sightingdata$yrs[i]))]<-sightingdata$sights[i]}
	d1<-which(res$sights>0)
	d2<-d1[length(d1)]
	test.dat<-res[(d2+1):length(res[,1]),]
	test.dat$chance<-0
	for(g in 1:length(test.dat$yrs)){
		#g=1
		dat<-res[res$yrs<=test.dat$yrs[g],]
		test.dat$chance[g]<-Solow1993.eq2.fun(dat)
	}
	cut<-subset(test.dat, test.dat$chance<=alpha)
	fin<-data.frame(Estimate=cut$yrs[1])
		if(data.out==T){
		return(test.dat[,c(1,3)])}
	else{return(fin)}	
}
