`print.scad.svm` <-
function(x,...){
	cat("\nSCAD-SVM\n")
	cat("\nBias = ", x$b)
	#cat("\nSelected Variables= ", colnames(x$x)[x$xind])#
	cat("\nSelected Variables= ", names(x$w))
	cat("\nCoefficients:\n  ")
	print(x$w)
	cat("\n")
}

