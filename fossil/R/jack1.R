`jack1` <-
function(x,taxa.row=TRUE,abund=TRUE) {
  if (taxa.row==FALSE) x<-t(x) 
  if (ncol(as.matrix(x))==1|nrow(as.matrix(x))==1) x <- x 
  else if (abund==TRUE) {x <- rowSums(x) 
    m <- sum(x)}
  else {m <-length(x[1,])
    x1 <- numeric(length(x[,1]))
    for (i in 1:length(x[,1])) x1[i] <- length(x[i,][x[i,]>0])
    x<-x1}
q1<- length(x[x==1])
m <- sum(x)
so <- length(x[x>0])
foj <- so+q1*((m-1)/m)
return(foj)
}

