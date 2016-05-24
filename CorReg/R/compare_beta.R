# ' compare signs of the coefficients in two vectors
# ' @export
# ' @param truebeta first vector
# ' @param hatbeta second vector
# ' @param taux boolean. Computes the ratio of each selection statistic or not.

compare_beta<-function(truebeta=truebeta,hatbeta=hatbeta,taux=FALSE){
   resloc=compare_zero(trueA=truebeta,Aalgo=hatbeta,taux=FALSE)
   resloc=c(resloc,compare_sign(trueA=truebeta,Aalgo=hatbeta))
   resloc$biais=sum((truebeta-hatbeta)^2)
   return(resloc)  
}


