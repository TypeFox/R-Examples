plot.crackRparallel <-
function(x, thresh=NA, sfpof.int=FALSE, ...)
{

    obj <- x[[1]]
    comb <- x[[2]]
    rm(x)

    if( length(obj) <= 1 ) stop("need two or more sets of results to plot crackRparallel object")

    if( sfpof.int )
    {
        comb$sfpof <- calcSfpofFromPofInt(comb$pof.int)
        for( iii in 1:length(obj) )
            obj[[iii]]$sfpof <- calcSfpofFromPofInt(obj[[iii]]$pof.int)
    }

    ## if sfpof has zero values, we need to provide a small positive number for plotting
    for( iii in 1:length(obj) )
        if(min(obj[[iii]]$sfpof$sfpof)<=0)
            obj[[iii]]$sfpof$sfpof[obj$sfpof$sfpof <= 0] <- min(obj[[iii]]$sfpof$sfpof[obj[[iii]]$sfpof$sfpof>0])

    plot.crackRresults(comb, col=4, boot=FALSE, ...)
    for(iii in 1:length(obj))
        lines(x=obj[[iii]]$sfpof$flight, y=obj[[iii]]$sfpof$sfpof, ...)

    mtext(paste("Number of parallel runs =", length(obj)),
          side=3, adj=1, cex=0.7, line=0.5)

    mtext("Combined results in blue",
          side=3, adj=1, cex=0.7, col=4, line=1.3)
    
    if( !is.na(thresh) )
        {
            abline(h=thresh, lty=3, lwd=2, col="red")
            mtext(paste("SFPOF threshold =", thresh),
                  side=3, adj=1, cex=0.7, line=1.3, col="red")
        }
}
