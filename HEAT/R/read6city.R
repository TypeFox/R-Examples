read6city <-
function(data, code){

citydata=subset(data, data[,1]==code)
ncc=dim(data)[2]
nrr=dim(data)[1]

cbind(citydata)
}
