predict.gppr<-function(object,newdata,type="link"){
	pred<-predict(object$ppr,newdata)
	if(type=="response"){
		pred<-object$family$linkinv(pred)
	}
	pred
}
