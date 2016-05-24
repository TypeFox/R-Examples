varwiselogLik.regression <-
function(xx, D, glm.family) {
	n<-nrow(xx)
	fit<-vector("list", n)
	for(i in 1:n) {
		# glm des vollen/reduzierten Modells
		fit[[i]]<-glm.fit(x=D, y=t(xx)[, i], family=glm.family)
	}
	return(fit)
}
