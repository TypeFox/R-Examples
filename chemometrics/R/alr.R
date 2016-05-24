"alr" <-
function(X,divisorvar){
# additive logratio transformation
# INPUT: X...data matrix, divisorvar...number of ratioing variable
X.alr <- log10(X/X[,divisorvar])
return(X.alr[,-divisorvar])
}

