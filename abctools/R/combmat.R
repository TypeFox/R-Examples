combmat <-function (n, limit = NULL) {


if (is.null(limit)) {
	limit <- n
}
    
limit <- min(limit, n)
m <- NULL
        
for (i in 1:limit) {
	tmp <- fillcomb(combn(n, i))
        m <- rbind(m, tmp)
}
    
colnames(m)<-paste("C",1:n,sep="")

m
}

