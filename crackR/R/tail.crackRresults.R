tail.crackRresults <-
function(x, sfpof.int=FALSE, ...)
{
    ## for sfpof.int==TRUE, use pof.int to approx SFPOF in each interval
    if( sfpof.int ) x$sfpof <- calcSfpofFromPofInt(x$pof.int)

    x$sfpof   <- tail(x$sfpof)
    x$pcd     <- tail(x$pcd)

    try(x$pof.int <- tail(x$pof.int), silent=TRUE)

    try(x$MC.results$fail.count <- tail(x$MC.results$fail.count), silent=TRUE)
    try(x$MC.results$repair.count <- tail(x$MC.results$repair.count), silent=TRUE)
    try(x$MC.results$all.results <- tail(x$MC.results$all.results), silent=TRUE)
    
    print(x, ...)
}
