"P.chess03" <-
function (filename="chess03.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=6.5, width=4)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	chessboard(seed=26)

        if(length(filename)) dev.off()
}

