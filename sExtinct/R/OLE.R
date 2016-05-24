OLE <-
function(sightingdata, alpha){

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
      
 	OLE.fun(sightingdata, alpha)
		}
