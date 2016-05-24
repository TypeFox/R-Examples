model2index <-
function (model,dig) { 			## need to tell it number of columns in maximal design matrix
res<-t(sapply(X=model,FUN=hex2bin))			## uses hex2bin function from BMS

s<-dim(res)[2]
res[,(s-dig+1):s]}
