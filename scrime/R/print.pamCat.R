`print.pamCat` <-
function(x,digits=3,...){
	cat("Prediction Analysis of Categorical Data\n\n")
	mat.theta<-x$mat.theta
	mat.theta[,3]<-round(mat.theta[,3],digits)
	print(mat.theta)
}

