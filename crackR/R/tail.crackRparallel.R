tail.crackRparallel <-
function(x, sfpof.int=FALSE, ...)
{
    for(iii in 1:length(x))
        {
            cat(paste("\n##\ntail() of run ", iii, "\n##\n\n", sep=""))
            tail(x[[iii]], sfpof.int, ...)
        }
}
