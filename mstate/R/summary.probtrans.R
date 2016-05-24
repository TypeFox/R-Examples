summary.probtrans <- function(object,from,complete=FALSE,variance=TRUE,...)
{
    if (!inherits(object, "probtrans"))
        stop("'object' must be a 'probtrans' object")
    cat("An object of class 'probtrans'\n")
    obj1 <- object[[1]]
    trans <- object$trans
    S <- nrow(trans)
    tt <- unique(obj1$time) # the time points
    nt <- length(tt)
    if (missing(from)) from <- 1:S
    startingStates <- intersect(from,1:S)
    if (nt<=12 | complete) {
        for (s in startingStates) {
            objs <- object[[s]]
            cat("\nPrediction from state",s,":\n\n")
            if (variance) print(objs,...)
            else print(objs[,1:(S+1)],...)
        }
    }
    else {
        for (s in startingStates) {
            objs <- object[[s]]
            cat("\nPrediction from state",s,"(head and tail):\n\n")
            if (variance) {
                print(head(objs),...)
                cat("\n...\n\n")
                print(tail(objs),...)
            }
            else {
                print(head(objs[,1:(S+1)]),...)
                cat("\n...\n\n")
                print(tail(objs[,1:(S+1)]),...)
            }
        }
    }
    return(invisible())
}
