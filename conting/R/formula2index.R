formula2index <-
function(big.X,formula,data){

dX<-model.matrix(formula,data=data,contrasts =attributes(big.X)$contrasts)## design matrix given by formula

indo<-c()
for(i in 1:dim(big.X)[2]){
indo[i]<-ifelse(any((dimnames(big.X)[[2]])[i]==(dimnames(dX)[[2]])),1,0)}	## which columns of maximal design matrix 											## are in this new matrix
indo}
