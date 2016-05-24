phi <-
function(X1,Y1,X2,Y2,Dmat,delta){

if (delta!=0){

phi1<-0.5*(((abs(X1-Y1))**delta)%*%Dmat%*%((abs(X1-Y1))**delta)+((abs(X2-Y2))**delta)%*%Dmat%*%((abs(X2-Y2))**delta))
phi2<-0.5*(((abs(X1-Y2))**delta)%*%Dmat%*%((abs(X1-Y2))**delta)+((abs(X2-Y1))**delta)%*%Dmat%*%((abs(X2-Y1))**delta))

return(c(phi1,phi2))
}
if (delta==0){
phi1<-0.5*(((abs(X1-Y1))!=0)%*%Dmat%*%((abs(X1-Y1))!=0)+((abs(X2-Y2))!=0)%*%Dmat%*%((abs(X2-Y2))!=0))
phi2<-0.5*(((abs(X1-Y2))!=0)%*%Dmat%*%((abs(X1-Y2))!=0)+((abs(X2-Y1))!=0)%*%Dmat%*%((abs(X2-Y1))!=0))

return(c(phi1,phi2))


}
}

