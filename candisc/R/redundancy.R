# calculate redundancy indices for canonical correlation analysis

redundancy <- function(object, ...) {
    if (!inherits(object, "cancor")) 
        stop("Not a cancor object")
    cancor <- object$cancor
    Xstruc <- object$structure$X.xscores
    Ystruc <- object$structure$Y.yscores
    
    # for each canonical variate, fraction of total X, Y variance associated
    Xcan.vad <- apply(Xstruc^2, 2, mean, na.rm = TRUE)
    Ycan.vad <- apply(Ystruc^2, 2, mean, na.rm = TRUE)
    
    # canonical redundancies for X, Y variables (total fraction of X variance accounted for by Y variables through canonical
    # variables, and vice-versa)
    Xcan.redun <- Xcan.vad * cancor^2
    Ycan.redun <- Ycan.vad * cancor^2
    
    result <- list(Xcan.redun=Xcan.redun, 
	               Ycan.redun=Ycan.redun, 
	               X.redun=sum(Xcan.redun), 
	               Y.redun=sum(Ycan.redun),
	               set.names=object$names$set.names)
    class(result) <- "cancor.redundancy"
    # invisible(result)
    result
} 

print.cancor.redundancy <- function(x, digits=max(getOption("digits") - 3, 3), ...) {
	Xname <- x$set.names[1]
	Yname <- x$set.names[2]
	cat(paste("\nRedundancies for the", Xname, "variables & total X canonical redundancy\n\n"))
	Xredun <- c(x$Xcan.redun, "total X|Y"=x$X.redun)
	print(Xredun, digits=digits)
	
	cat(paste("\nRedundancies for the", Yname, "variables & total Y canonical redundancy\n\n"))
	Yredun <- c(x$Ycan.redun, "total Y|X"=x$Y.redun)
	print(Yredun, digits=digits)
	
}
