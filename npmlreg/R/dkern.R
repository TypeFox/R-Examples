"dkern" <-
function(x,y,k,lambda){
    ifelse(y==x, lambda, (1-lambda)/(k-1))
}

