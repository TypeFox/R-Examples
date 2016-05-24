print.crackRparallel <-
function(x, sfpof.int=FALSE, ...)
{
    for(iii in 1:length(x))
        {
            cat(paste("\nrun ", iii, "\n##\n\n", sep=""))
            print(x[[iii]], sfpof.int, ...)
        }
}
