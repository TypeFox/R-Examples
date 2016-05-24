"P.sticks03" <-
function (filename="sticks03.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=3, width=5)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	sticks(seed=9)

        if(length(filename)) dev.off()
}

