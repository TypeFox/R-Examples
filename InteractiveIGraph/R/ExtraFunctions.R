`%.%` <- function (x, y) paste(x, y, sep = "")

d2file <- function(file="Rplot%d",height=7.5,width=7.5,which=dev.cur(), ...) {
	if(which!=dev.cur()) {
		cur <- dev.cur()
		dev.set(which)
	}
	else cur <- NULL

	dev.copy(device=pdf,file=file,height=height,width=width,
		onefile=TRUE,paper="special",...)	
	dev.off()
	if(!is.null(cur))dev.set(cur)
	else dev.set(which)
}