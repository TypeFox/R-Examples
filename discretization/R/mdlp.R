mdlp <-
function(data){
    p <- length(data[1,])-1
    y <- data[,(p+1)]
    xd <- data
    cutp <- list()
    for (i in 1:p){
        x <- data[,i]
        cuts1 <- cutPoints(x,y)
        cuts <- c(min(x),cuts1,max(x))
        cutp[[i]] <- cuts1
        if(length(cutp[[i]])==0)cutp[[i]] <- "All"
        xd[,i] <- as.integer(cut(x,cuts,include.lowest = TRUE))
    }
    return (list(cutp=cutp,Disc.data=xd))
}
