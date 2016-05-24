"P.chess" <-
function (filename="chess.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=6.5, width=4)
                par(mar=c(1,0,1,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	chessboard(seed=13)

        if(length(filename)) dev.off()
}

