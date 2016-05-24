"P.polytile03" <-
function (filename="polytile03.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=3, width=5)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	polytile(seed=7)

        if(length(filename)) dev.off()
}

