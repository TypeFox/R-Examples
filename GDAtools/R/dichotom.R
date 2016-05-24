dichotom <- function(data,out='numeric') {
  if(!is.data.frame(data)) data <- data.frame(data)
  res <- matrix(nrow=nrow(data),ncol=length(levels(data[,1])))
  for(i in 1:ncol(data)) {
    if(is.factor(data[,i])==FALSE) data[,i] <- factor(data[,i])
    nlevels <- length(levels(data[,i]))
    temp <- matrix(nrow=nrow(data),ncol=nlevels)
    for(j in 1:nlevels) temp[,j] <- ifelse(data[,i]==levels(data[,i])[j],1,0)
    colnames(temp) <- paste(names(data)[i],levels(data[,i]),sep=".")
    if(i==1) res <- temp else res <- cbind(res,temp)
    }
  res <- as.data.frame(res)
  if(out=='factor') for(i in 1:ncol(res)) res[,i] <- as.factor(res[,i])
  res
}