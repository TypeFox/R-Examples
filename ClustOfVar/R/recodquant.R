
recodquant <-
  function(X)
  {
    X <- as.matrix(X)
    Xcod <- apply(X,2,missing.mean)
    red <- sqrt((nrow(X)-1)/nrow(X))
    sd.Xcod <- apply(Xcod,2,sd)*red
    mean.Xcod <- apply(Xcod,2,mean)
    Z<- scale(Xcod,scale=sd.Xcod) 
    apply(Z,1,function(x) sum(is.na(x))) 
    if (sum(is.na(Z))!= 0) stop("There are columns in X.quanti where all the values are identical")
    return(list(Z=Z,g=mean.Xcod,s=sd.Xcod,Xcod=Xcod))
  }