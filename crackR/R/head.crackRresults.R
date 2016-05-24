head.crackRresults <-
function(x, sfpof.int=FALSE, ...)
{
    ## for sfpof.int==TRUE, use pof.int to approx SFPOF in each interval
    if( sfpof.int ) x$sfpof <- calcSfpofFromPofInt(x$pof.int)

    x$sfpof   <- head(x$sfpof)
    x$pcd     <- head(x$pcd)

    try(x$pof.int <- head(x$pof.int), silent=TRUE)

    try(x$MC.results$fail.count <- head(x$MC.results$fail.count), silent=TRUE)
    try(x$MC.results$repair.count <- head(x$MC.results$repair.count), silent=TRUE)
    try(x$MC.results$all.results <- head(x$MC.results$all.results), silent=TRUE)
    
    print(x, ...)
}
