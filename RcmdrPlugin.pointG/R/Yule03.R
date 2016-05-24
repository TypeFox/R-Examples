Yule03<-
function (X,YY,levX,nameYY,tri=0) 
{

    L <- ncol(YY)
    if (is.factor(YY[,1])) {
        T1 <- Yule01(X,YY[,1],levX,nameYY[1])
    }
    else {
        T1 <- Yule02(X,YY[,1],levX,nameYY[1])
    }
    if (L > 1) {
        for (l in 2:L) {
            if (is.factor(YY[,l])) {
                T2 <- Yule01(X,YY[,l],levX,nameYY[l])
            }
            else {
                T2 <- Yule02(X,YY[,l],levX,nameYY[l])
            }
            T1 <- rbind(T1, T2)
        }
    }
T1[,2]<-round(T1[,2],2)
T1[,3]<-round(T1[,3],2)

    if(tri==0){
return(T1)}
else{
return(T1[order(T1[,2]),])}
}