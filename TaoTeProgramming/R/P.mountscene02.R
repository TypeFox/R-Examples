"P.mountscene02" <-
function (filename="mountscene02.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=6.5, width=4)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	mountscene(seed=10, df=c(seq(4, 40, length=6), 40, 30, 20, 20), 
		levels=10)

        if(length(filename)) dev.off()
}

