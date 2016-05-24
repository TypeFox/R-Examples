## compute AUC
AUC <- function(x, y, group, switchAUC = TRUE){
    if(missing(y) && missing(group)){
        stop("'y' or 'group' must be specified!")
    }
    if(missing(y) && !missing(group)){
        group <- factor(group)
        if(nlevels(group) != 2)
            stop("'group' must have exactly two different group labels!")
        group <- as.integer(group)-1
        y <- x[group == 0]
        x <- x[group == 1]
    }
    x <- sort(x)
    y <- sort(y)
    nx <- length(x)
    ny <- length(y)
    AUC <- (nx*ny + nx*(nx+1)/2 - sum(rank(c(x,y))[1:nx]))/(nx*ny)
    if(AUC < 0.5 & switchAUC){
	dig <- getOption("digits")
	warning("The computed AUC value ", round(AUC, dig), " will be replaced by ", 
		round(1-AUC, dig),  " which can be achieved be interchanging the sample labels!")
        return(1-AUC)
    }else
        return(AUC)
}
