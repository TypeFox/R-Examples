`bootstrap` <-
function(x,taxa.row=TRUE,abund=TRUE,samples=NA) {
  if (taxa.row==FALSE) x<-t(x) 
  if (ncol(as.matrix(x))==1|nrow(as.matrix(x))==1) {x <- x 
  m<-samples}
  else if (abund==TRUE) {x <- rowSums(x) 
    m <- sum(x)}
  else {m <-length(x[1,])
    x1 <- numeric(length(x[,1]))
    for (i in 1:length(x[,1])) x1[i] <- length(x[i,][x[i,]>0])
    x<-x1}
  so <- length(x[x>0])
  spk<-0
  for (i in 1:so) spk<-spk+((1-x[i]/m)^m)
  sboot <- so+spk
  attr(sboot,"method") <- "Bootstrap"
  return(sboot)
}

