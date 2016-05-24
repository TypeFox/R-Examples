# IDENTIFY CONSTRAINTS FOR A GIVEN CRITERION
which.at.boundary <- function(object, criterion = 1e-6){

if(class(object)[1]!="lexpit"&class(object)[1]!="blm")
		stop("Object must be an instance of a blm or lexpit model.")
	
if(any(predict(object)<=criterion|predict(object)>=1-criterion)){
	which <- which(predict(object)<=criterion|predict(object)>=1-criterion)
	if(class(object)[1]=="lexpit")
		cbind(model.matrix(object@formula.linear,object@data)[,-1],
				 model.matrix(object@formula.expit,object@data)[,-1])[which,]
	else
		cbind(model.matrix(object@formula,object@data)[,-1])[which,]
}
else{
	cat("No boundary constraints using the given criterion.\n\n")
	invisible(NA)
}
		
}