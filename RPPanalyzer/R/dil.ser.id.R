`dil.ser.id` <-
function(x,sample.id=c("sample_name","sample_treat")){
               
        mmnt.lines <- which ( x[[4]][,"sample_type"]=="measurement")
	     temp <- c(NULL)
   		for (i in sample.id){
   
   			temp <- paste(temp,as.character(x[[4]][mmnt.lines,i]),sep="")
   
   			}
   		identifier <- unique(temp)
   
   				dilution.steps <- (length(temp)/length(identifier))
   				count <- length(identifier)
   
   			print(paste("found",count,"individual samples, spottet in",dilution.steps
   					,"step dilution series"))
   	 
        
   	 return(identifier)
}

