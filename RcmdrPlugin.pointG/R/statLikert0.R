statLikert0<-function (X,tri=0,rr=2) 
{
     
        X <- as.matrix(X)
        Mean <- apply(X, 2, mean, na.rm = TRUE)
        Min <- min(X, na.rm = TRUE)
        Max <- max(X, na.rm = TRUE)
        if (tri==0) {
 return(round(Mean,rr))
    }
    else {
return(round(sort(Mean),rr))
           }
}