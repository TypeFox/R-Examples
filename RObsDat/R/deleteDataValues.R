deleteDataValues <- function(ID=NULL,reason=NULL){
  
	the.class <- class(ID)
	the.id=c()
	
	if(any(the.class=="inherited_stfdf")){
		object = ID
		# convert data from data.frame to vector
		for(column in ncol(object@ValueIDs@data)){
			the.id = c(the.id, object@ValueIDs@data[,column])
	        }
	}  else if(any(the.class=="STFDF")){
		object = ID
		# convert data from data.frame to vector
		for(column in ncol(object@data)){
			the.id = c(the.id, object@data[,column])
	        }
	}  else if(any(the.class=="xts")){
		the.id = coredata(ID)
	} else if(the.class=="integer"){
		the.id = ID      
	} else {
		#warning("If only one time point is selected you have to use '@ValueIDs@data' in addition.")
		warning("Not sure this should happen")
	}
	
	
	if(length(the.id)>0){
		IarchiveDataValues(getOption("odm.handler"),ValueID=the.id, reason)
		IdeleteDataValues(getOption("odm.handler"),ValueID=the.id)
	}

}
