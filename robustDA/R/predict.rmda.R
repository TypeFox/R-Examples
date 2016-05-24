predict.rmda <-
function(object,X,...){
# Predict classes for new data
	prms = object$prms
	R = object$R
	P = estep(prms$modelName, X, prms$parameters)$z
	Pfin = P %*% t(R)
	cls = max.col(Pfin)
	res = list(cls=cls,P=Pfin)
}
