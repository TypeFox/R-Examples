disc.Topdown <-
function(data, method=1){
    type <- c("CAIM", "CACC", "ameva")
    meth <- type[method]
    cutList <- topdown(data,method)
    p <- length(data[1,])-1
    xd <- data
    for (i in 1:p){
        cuts <- cutList[[i]]
        xd[,i] <- as.data.frame(as.integer
                    (findInterval(data[,i],cuts,rightmost.closed=TRUE)))
    }
    list(cutp=cutList,Disc.data=xd)
}
