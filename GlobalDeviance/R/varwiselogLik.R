varwiselogLik <-
function(xx, D, glm.family) {
	n<-nrow(xx)
	fit<-vector("list", n)
	fit<-varwiselogLik.regression(xx=xx, D=D, glm.family=glm.family)
	# Extraktion log-Likelihood
	res<-matrix(NA, nrow=n, ncol=1)
	res[, 1]<-sapply(fit, function(x) {  deviance(x) })
	return(res)
}
