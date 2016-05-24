`AcfPlot` <-
function(g, LagZeroQ=TRUE, ylab=NULL, main=NULL, ...){
if (LagZeroQ) {
    plot( 0:length(g), c(1,g), type="h", ylim=c(-1,1), xlab="lag", ylab=ylab, main=main)
    lines( c(0, length(g)), c(0,0), col="magenta")
    }
else {
    plot( 1:length(g), g, type="h", ylim=c(-1,1), xlab="lag", ylab=ylab, main=main)
    lines( c(0, length(g)), c(0,0), col="magenta")
    }
}

