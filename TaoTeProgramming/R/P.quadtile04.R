"P.quadtile04" <-
function (filename="quadtile04.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=3, width=5)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	quadtilebalance(c(3,5), seed=32)

        if(length(filename)) dev.off()
}

