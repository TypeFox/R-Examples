Dint <- function(x,y,range,restrictNonMissing=FALSE){
  if (is.null(range)) range=c(x[1],x[length(x)])
  ##   integrate a step function f with
  ##   values y=f(x) between range[1] and range[2]
  start <- max(range[1],min(x))
  Stop <- min(range[2],max(x))
  if ((Stop-start)<=0)
    return(0)
  else{
    Y=y[x>=start & x<Stop]
    X=x[x>=start & x<Stop]
    if (restrictNonMissing){
      X=X[!is.na(Y)]
      Y=Y[!is.na(Y)]
    }
    else if (any(is.na(Y))|| any(is.na(X))){
      return(NA)
    }
    return(1/(Stop-start) * sum(Y*diff(c(X,Stop))))
  }
}
