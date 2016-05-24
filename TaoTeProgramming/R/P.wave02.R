"P.wave02" <-
function (filename="wave02.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=6.5, width=4)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	waves02(lwd=10, seed=7)

        if(length(filename)) dev.off()
}

