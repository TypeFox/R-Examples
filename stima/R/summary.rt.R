summary.rt <-
function(object,digits=3,...){
	object$goffull<-round(object$goffull,digits=digits)
	if(object$goffull[1,5]==0){object$goffull<-object$goffull[,-c(5:6)]}

	if(is.null(object$gofsel)){
	return(object$gof)}
	else{
		if(object$gofsel[3]==0){object$gofsel<-object$gofsel[-c(3:4)]}
		gofsel<-round(object$gofsel,digits=digits)
	   list(full=object$goffull,selected=gofsel)}
}
