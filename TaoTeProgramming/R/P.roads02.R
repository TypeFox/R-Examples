"P.roads02" <-
function (filename="roads02.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=6.5, width=4)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	roads(100, seed=3)

        if(length(filename)) dev.off()
}

