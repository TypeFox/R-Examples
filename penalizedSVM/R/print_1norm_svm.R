`print.1norm.svm` <-
function(x,...){
	cat("\n 1norm-SVM\n")
	cat("\nBias = ", x$b)
	cat("\nSelected Variables= ", names(x$w))
	cat("\nCoefficients:\n  ")
	print(x$w)
}

