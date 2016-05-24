#'@method print plant3dlist
#'@S3method print plant3dlist
print.plant3dlist <- function(x, ...){
	cat("Yplant - object of class 'plant3dlist'.\n\n")
	cat(paste(c(rep("-",30),"\n"),collapse=""))
	cat("Number of plants included :", attributes(x)$nplants, "\n")
	cat("File names:\n\n")
	print(data.frame(Pfile=attributes(x)$pfiles, Lfile=attributes(x)$lfiles))
}
