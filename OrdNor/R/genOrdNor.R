genOrdNor <-
function(n,plist, cmat.star, mean.vec, sd.vec, no.ord, no.norm) {

if ( length(mean.vec) != no.norm ){
stop("Dimension of the mean vector does not match the number of normal variables!\n")
} 

if ( length(sd.vec) != no.norm ){
stop("Dimension of the sd vector does not match the number of normal variables!\n")
} 

if (min(sd.vec)<0){
stop("Standard deviation values cannot be negative!\n")
} 
if (no.norm==0){
YY = ordsample(n, plist, cmat.star, cormat="continuous")
}

if (no.ord==0){
YY = rmvnorm(n, rep(0,ncol(cmat.star) ), cmat.star)
YY = t( t(YY)*sd.vec+mean.vec )
}

if (no.norm>0 & no.ord>0) {
XX=rmvnorm(n, rep(0,ncol(cmat.star) ), cmat.star)
YY=NULL;
for( i in 1:length(plist)  ) {
OO = ordinalize(plist[[i]], XX[,i])
YY=cbind(YY, OO) 
rm(OO)
}

YY = cbind(YY, XX[, (length(plist)+1) : ncol(cmat.star) ] )
YY[,(length(plist)+1) : ncol(cmat.star)] = t( t(YY[,(length(plist)+1) : ncol(cmat.star)])*sd.vec+mean.vec )
rm(XX)
}
colnames(YY)<-NULL
return(YY)
}
