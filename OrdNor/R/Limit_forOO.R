Limit_forOO <-
function(pvec1, pvec2){

validate.plist(list(pvec1),1)
validate.plist(list(pvec2),1)

X=rnorm(100000,0,1)
Y=rnorm(100000,0,1)

XORD = ordinalize(pvec1, X)
YORD = ordinalize(pvec2, Y)

max = cor(XORD[order(XORD)],YORD[order(YORD)])
min = cor(XORD[order(XORD,decreasing=TRUE)],YORD[order(YORD)])

rm(X,Y)
return(c(min,max))
}
