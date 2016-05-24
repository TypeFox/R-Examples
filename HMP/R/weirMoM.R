weirMoM <-
function(data, MoM, se=FALSE){
K <- ncol(data)
J <- nrow(data)
rowSumsdata <- rowSums(data) + 0.000001
colSumsData <- colSums(data)

MSP <- (J-1)^(-1) * sum(rowSums((data/rowSumsdata - matrix(rep(MoM, J), J, K, byrow=TRUE))^2) * rowSumsdata)
MSG <- (sum(colSumsData)-J)^(-1) * sum(rowSums(data/rowSumsdata * (1-data/rowSumsdata)) * rowSumsdata)
nc <- 1/(J-1) * (sum(rowSumsdata)-sum(rowSumsdata^2)/sum(rowSumsdata))
MoM.wh <- (MSP-MSG)/(MSP+(nc-1)*MSG)

if(se){
std.er <- sqrt(2 * (1-MoM.wh)^2/(J-1) * ((1+(nc-1) * MoM.wh)/nc)^2)
return(list(theta=MoM.wh, se=std.er))
}else{
return(MoM.wh)
}
}
