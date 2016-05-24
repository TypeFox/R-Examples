"P.backcolor" <-
function (filename="backcolor.png") 
{
        if(length(filename)) {
                png(filename=filename, height=2500, width=1838)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	red <- c(250, 168, 176, 184, 192, 200, 208)
	green <- c(255, 168, 176, 184, 192, 200, 208)
	blue <- c(250, 208, 216, 224, 232, 240, 248)
	mountscene(seed=11, df=c(seq(4, 40, length=4), 30, 20, 20), levels=7,
		box=FALSE, color=rgb(red, green, blue, maxColorValue=255))
	points(runif(500), runif(500), pch=".", col="white", cex=2)
	points(runif(100), runif(100), pch=".", col="white", cex=5)

        if(length(filename)) dev.off()
}

