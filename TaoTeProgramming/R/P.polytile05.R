"P.polytile05" <-
function (filename="polytile05.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=5.3, width=4)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	polytile(seed=14)

        if(length(filename)) dev.off()
}

