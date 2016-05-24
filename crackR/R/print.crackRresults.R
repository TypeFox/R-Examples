print.crackRresults <-
function(x, sfpof.int=FALSE, ...)
{
    ## for sfpof.int==TRUE, use pof.int to approx SFPOF in each interval
    if( sfpof.int ) x$sfpof <- calcSfpofFromPofInt(x$pof.int)
    class(x) <- "list"
    print(x)
}
