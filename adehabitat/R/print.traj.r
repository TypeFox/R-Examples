"print.traj" <- function(x, ...)
{
    if (!inherits(x, "traj"))
        stop("x should be an object of class traj")
    levani<-levels(x$id)
    u<-split(x$burst, x$id)
    cat("******** Data frame of class traj *********\n\n")
    for (i in 1:length(u)) {
        cat("Animal ",names(u)[i],":   ",
            nlevels(factor(u[[i]])), " circuits")
        cat(" (",length(u[[i]])," relocations)\n", sep="")
    }
    cat("\nVariables measured for each relocation:\n\n")
    print(names(x), quote=FALSE, ...)
}

