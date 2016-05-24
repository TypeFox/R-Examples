"P.frontcolor" <-
function (filename="frontcolor.png") 
{
        if(length(filename)) {
                png(filename=filename, height=2500, width=1838)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	red <- c(172, 160, 148, 136, 124, 112, 0)
	green <- c(172, 160, 148, 136, 124, 112, 140)
	blue <- c(252, 240, 228, 216, 204, 192, 0)
	mountscene(seed=11, df=c(seq(4, 40, length=4), 30, 20, 20), levels=7,
		box=FALSE, color=rev(rgb(red, green, blue, maxColorValue=255)))

        if(length(filename)) dev.off()
}

