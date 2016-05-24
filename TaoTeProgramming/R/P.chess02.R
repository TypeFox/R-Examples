"P.chess02" <-
function (filename="chess02.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=7, width=5)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	chessboard(seed=39)

        if(length(filename)) dev.off()
}

