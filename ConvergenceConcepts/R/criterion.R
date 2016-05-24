# Comnputation of the criterion from the sample paths
#####################################################

criterion <- function(data,epsilon=0.05,mode="p",r=2) {


  M <- nrow(data)

  if (mode=="p") crit <- colSums(abs(data)>epsilon)/M 

  if (mode=="as") crit <- rev(colSums(t(apply(t(apply(abs(data)>epsilon,FUN=rev,MARGIN=1)),
	FUN=cumsum,MARGIN=1))>0)/M) 
  
  if (mode=="r") {data <- abs(data); data <- data^r; crit <- colMeans(data)} 
  
  return(list(crit=crit))

}

