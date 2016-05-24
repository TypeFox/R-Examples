Limit_forON <-
function(pvec1){
validate.plist(list(pvec1),1)

X=rnorm(100000,0,1)
Y=rnorm(100000,0,1)

XORD = ordinalize(pvec1, X)

max = cor(XORD[order(XORD)],Y[order(Y)])
min = cor(XORD[order(XORD,decreasing=TRUE)],Y[order(Y)])

rm(X,Y)
return(c(min,max))
}
