"P.polytile04" <-
function (filename="polytile04.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=4, width=4)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	polytile(seed=10)

        if(length(filename)) dev.off()
}

