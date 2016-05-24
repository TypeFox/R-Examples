#'@method print projectedplant3d
#'@S3method print projectedplant3d
#'@rdname projectplant
print.projectedplant3d <- function(x, ...){
	cat("Yplant - object of class 'projectedplant3d'.\n\n")
	cat(paste(c(rep("-",30),"\n"),collapse=""))
	cat("Plant has",length(x$leaves),"leaves.\n")
	cat("Projection from azimuth =",x$viewangle$azimuth,"and altitude =",x$viewangle$altitude,"\n")
}
