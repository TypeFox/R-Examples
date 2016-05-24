"P.tritile02" <-
function (filename="tritile02.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=4, width=4)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	tritile(seed=8)

        if(length(filename)) dev.off()
}

