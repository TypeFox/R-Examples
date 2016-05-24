predict.sublasso<-function(object, xpred, type, s, ...){
 if (missing(type)) {type="class"}
 if (missing(s))    {s=object$lambda}
 predy=predict(object$fit,t(xpred),s=s,type=type)
 return(list(predy=predy))
}


