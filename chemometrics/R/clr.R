"clr" <-
function(X){
# centered logratio transformation
Xgeom <- exp(1)^apply(log(X),1,mean)
X.clr <- log(X/Xgeom)
return(X.clr)
}

