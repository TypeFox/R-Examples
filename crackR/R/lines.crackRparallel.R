lines.crackRparallel <-
function(x, sfpof.int=FALSE, ...)
{
    for(iii in 1:length(x))
    {
        if( sfpof.int ) x[[iii]]$sfpof <- calcSfpofFromPofInt(x[[iii]]$pof.int)
        lines(x[[iii]]$sfpof$flight, x[[iii]]$sfpof$sfpof, ...)
    }
}
