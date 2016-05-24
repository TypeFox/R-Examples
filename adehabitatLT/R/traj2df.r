".traj2df" <- function(x) {
    if (!inherits(x, "traj"))
        stop("x should be of class traj")
    class(x)<-"data.frame"
    row.names(x)<-as.character(1:nrow(x))
    return(x)
}

