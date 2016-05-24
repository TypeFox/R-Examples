plot.CD <-
function(x, sfpof.int=FALSE, ...)
{
    obj <- x
    rm(x)
    
    results     <- obj$results
    new.results <- obj$new.results

    if( sfpof.int )
        new.results$sfpof <- calcSfpofFromPofInt(new.results$pof.int)
    
    if( dim(results$sfpof)[1] > 0 )
        {
            plot(results, sfpof.int, ...)
            
            ## if new.results is present (meaning more than one run is in the crackR object),
            ##   add the sfpof line for this, but not PCD
            if( !is.null(new.results) )
                {
                    lines(x=new.results$sfpof$flight, y=new.results$sfpof$sfpof, col=3, lty=3)
                    if( obj$parameters$bootstrap.sfpof )
                        {
                            legend("topleft", c("Cum SFPOF Est", "New SFPOF Est", "Cum Quantiles"),
                                   lty=c(1,2,3), col=c(1,"forest green",4))
                        } else {
                            legend("topleft", c("Cum SFPOF Est", "New SFPOF Est"), lty=c(1,2), col=c(1,3))
                        }
                }
        }
}
