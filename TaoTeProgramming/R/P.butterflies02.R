"P.butterflies02" <-
function (filename="butterflies02.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=4.2, width=4)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	butterflies(seed=33)

        if(length(filename)) dev.off()
}

