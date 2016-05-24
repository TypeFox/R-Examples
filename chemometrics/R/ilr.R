"ilr" <-
function(X){
# isometric logratio transformation
X.ilr=matrix(NA,nrow=nrow(X),ncol=ncol(X)-1)
  for (i in 1:ncol(X.ilr)){
    X.ilr[,i]=sqrt((i)/(i+1))*log(((apply(as.matrix(X[,1:i]), 1, prod))^(1/i))/(X[,i+1]))
  }
return(X.ilr)
}

