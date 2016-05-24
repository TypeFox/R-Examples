modelPredict <- function(model,newdata){

	newdata 	   <- data.frame(newdata)
	names(newdata) <- names(model$data$X)
	
	if(is.null(model$type)){
		stop("Argument \'model\' should be a fitted model obtained from \'modelFit\'")
	} else type <- model$type
	
	if (type=="Linear" | type == "MARS" | type == "Additive" | type == "StepLinear"){
		y_pred <- predict(object=model$model,newdata=newdata)
	} else if (type=="PolyMARS"){
		y_pred <- polspline::predict.polymars(model$model,x=newdata)
	} else if (type=="Kriging"){
		y_pred <- predict(model$model,newdata,"UK")$mean
	} else stop("This method is not implemented yet.")
	return(as.numeric(y_pred))
}