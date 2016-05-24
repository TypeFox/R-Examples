getMahalStrata <-
function(col1, col2){
    col1 <- as.integer(col1)
    col2 <- as.integer(col2)
    dat <- data.frame(col1, col2)
    
    rownames(dat) <- 1:nrow(dat)
    g <- graph.data.frame(dat)
    links <- data.frame(col1=V(g)$name,group=clusters(g)$membership)
    return(merge(dat,links,by="col1",all.x=TRUE,sort=FALSE)$group)
}
