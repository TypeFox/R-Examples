head.crackRparallel <-
function(x, sfpof.int=FALSE, ...)
{
    for(iii in 1:length(x))
        {
            cat(paste("\n##\nhead() of run ", iii, "\n##\n\n", sep=""))
            head(x[[iii]], sfpof.int, ...)
        }
}
