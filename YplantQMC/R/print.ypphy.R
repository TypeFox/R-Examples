#'@method print ypphy
#'@S3method print ypphy
print.ypphy <- function(x,...){
	cat("Yplant - Leaf physiology object (class 'ypphy')\n\n")
    cat("Model :", x$leafmodel, "\n")
	cat("-------------------------------\n")
	cat("Parameters :\n\n")
	for(i in 1:length(x$leafpars)){
		p <- x$leafpars[i]
		cat(names(p),"=",p[[1]],"\n")
	}
}
