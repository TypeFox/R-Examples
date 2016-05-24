recodquant <-
  function(X)
  {
    X <- as.matrix(X)
    missing.mean <-
      function(C1){
        moy <- mean(C1,na.rm=T)
        ind <- which(is.na(C1)==T)
        if(length(ind)>=1){C1[ind]<-moy
        }
        return(C1)
      }
    Xcod <- apply(X,2,missing.mean)
    red <- sqrt((nrow(X)-1)/nrow(X))
    sd.Xcod <- apply(Xcod,2,sd)*red
    mean.Xcod <- apply(Xcod,2,mean)
    Z<- scale(Xcod,scale=sd.Xcod) 
    apply(Z,1,function(x) sum(is.na(x))) 
    if (sum(is.na(Z))!= 0) stop("There are columns in X.quanti where all the values are identical")
    return(list(Z=Z,g=mean.Xcod,s=sd.Xcod,Xcod=Xcod))
  }
