`TrenchMean` <-
function(r,z){
    g1<-colSums(TrenchInverse(toeplitz(r)))
    sum(g1*z)/sum(g1)
}

