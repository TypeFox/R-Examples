`InformationMatrixARp` <-
function(phi, lags){
p<-length(lags)
PMAX<-max(lags)
iar<-toeplitz(TacvfAR(phi, PMAX-1))
iar[lags,lags]
}

