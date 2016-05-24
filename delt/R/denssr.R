denssr<-function(volume,havlkm,n,method="loglik",mix=NULL){
#Calculates the log-likelihood 
#
if (havlkm==0) vastaus<-0
else if (method=="loglik")
                      vastaus<-havlkm*log(havlkm/(n*volume))
else if (method=="projec"){
   if (is.null(mix))  vastaus<-havlkm^2/(n*volume)
   else vastaus<-mix*(2-mix)*havlkm^2/(n^2*volume)
}
#
return(vastaus)
}



