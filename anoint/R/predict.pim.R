offset.pim <- function(object,data,index,beta,theta){
	
	X <- model.matrix(object@prognostic,data)
	if(all(X[,1]==1)) X <- X[,-1]
	S <- X%*%beta
	S[data[,index]==1] <- S[data[,index]==1]*theta
	
S
}

pim.predict <- function(object,...){
	
	X <- model.matrix(object@formula@formula@prognostic,object@formula@data)
	if(all(X[,1]==1)) X <- X[,-1]
	S <- X%*%object@coef[[2]]
	S[object@formula@data[,object@formula@trt.index]==1] <- 
		S[object@formula@data[,object@formula@trt.index]==1]*object@coef[[1]]
	
	if(length(object@coef[[1]])==1){
		S <- S+object@coef[[1]]*object@formula@data[,object@formula@trt.index]	
	}
	else{
		S <- S+object@coef[[1]][2]*object@formula@data[,object@formula@trt.index]
		S <- S+object@coef[[1]][1]
	}
S
}


setMethod("predict","pim",pim.predict)
