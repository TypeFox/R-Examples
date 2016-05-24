weightsit<-function(n,h,katka=4)
{
#normvakio<-(sqrt(2*pi)*h)^{-1}
resu<-matrix(0,n,1)
zumma<-0
for (i in 1:n){
    eta<-(n-i)
    if (eta/h>katka) tulos<-0 else tulos<-exp(-eta^2/(2*h^2))#*normvakio
    resu[i]<-tulos
    zumma<-zumma+tulos
}

resu<-resu/zumma

return(resu)
}



