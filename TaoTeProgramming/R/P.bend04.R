"P.bend04" <-
function (filename="bend04.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=6.5, width=4)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	bendplot(xdelta=100, ydelta=200, sd=.001, seed=5)

        if(length(filename)) dev.off()
}

