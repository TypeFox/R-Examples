missing.dist <-
function(data,...){	
      	if(is.null(dim(data))) n <- length(data) #data <- matrix(data)
	else n <- dim(data)[2]
	res <- NULL
	for(i in 1:n){
		id <- which(complete.cases(data[,i])==FALSE)
		if(length(id)>=1)res <- rbind(res,data.frame(i,id))
	}
	return(res)
  }
