ANGLE <- function(x,ref)
    {
    R <- sum(x*ref)/(sqrt(sum(x^2))*sqrt(sum(ref^2)))
    R <- ifelse(R < -1,-1,R)
    R <- ifelse(R > 1,1,R)
    acos(abs(R))
    }


anglesXYPlot <- function(x, angles="radiants", which=1:ncol(x$PPdir),...)
    {
    TYPE <- match.arg(angles, c("radiants", "degree"))
    
    X <- x$PPdir[,which]
    y <- X[,1]
    
    ANG <- apply(X,2,ANGLE,ref=y)
    
    if (TYPE=="degree") ANG <- ANG*180/pi
    
    ARGS <- list(...)
    #if (!hasArg(ylab)) ARGS <- c(list(ylab=paste("Angles in ",TYPE, sep="")),ARGS)
    #if (!hasArg(xlab)) ARGS <- c(list(xlab="Directions"),ARGS)
 
    ARGS <- c(list(ylab=paste("Angles in ",TYPE, sep="")),ARGS)
    ARGS <- c(list(xlab="Directions"),ARGS)
    
    do.call(xyplot, c(list(x=ANG~which),ARGS))
    }