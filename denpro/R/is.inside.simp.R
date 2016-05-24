is.inside.simp<-function(point,simple,rho)
{
d<-length(point)
dd<-1
sisalla2<-1
while ( (dd<=(d+1)) && (sisalla2==1) ){
      karki<-simple[dd,]
      dista2<-sum((point-karki)^2)
      if (dista2>rho^2){ 
            sisalla2<-0
      }
      dd<-dd+1
}

return(sisalla2)
}






