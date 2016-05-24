
writerandomdiffuse <- function(n, intensity=NULL){

	spher <- LeafAngle::angledist('spherical')
	zen <- LeafAngle::drawsample(spher,n)
	az <- runif(n,0,2*pi)
	
	if(is.null(intensity))intensity <- 1/n
	
	sky <- data.frame(zen=zen, az=az, intensity=intensity)
	write.table(sky, "randomdiffuse.dat",col.names=FALSE, row.names=FALSE)

}