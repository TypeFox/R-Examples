# ' Curve of the BIC for each possible p2 with a fixed Z and truncature of Z
# ' @param X matrix containing the dataset
# ' @param Z adjacency matrix (binary) describing the structure between the variables
# ' @param Bic_null_vect vector of the BIC for each variable. used when the variable is independant
# ' @param plot boolean to plot or not the curve
# ' @param star boolean to use BIC* (hierarchical uniform law on the structure)
# ' @param trunc number of sub-regression to keep (best R2). if NULL the min of BIC is kept
# ' @export
BicZcurve<-function(X=X,Z=Z,Bic_null_vect=Bic_null_vect,plot=T,star=F,trunc=NULL){
   p2=sum(colSums(Z)!=0)
   if(is.null(Bic_null_vect)){
      Bic_null_vect=density_estimation(X=X)$BIC_vect 
   }
   curve=sum(BicZ(X=X,Z=0*Z,Bic_null_vect=Bic_null_vect,star=star))
   
   if(p2>0){
      I2=which(colSums(Z)!=0)
      sigmavect=R2Z(Z=Z,X=X,crit="R2",adj=T)
      ordre=order(sigmavect[I2],decreasing=T)
      for (i in 1:p2){
         Zloc=Z;Zloc[,-I2[ordre[1:i]]]=0
         curve=c(curve,sum(BicZ(X=X,Z=Zloc,Bic_null_vect=Bic_null_vect,star=star)))
      }
      if(plot){
         plot(curve[-1])
         abline(h=curve[1])
      }
   }   
   quimin=which.min(curve)
   if(quimin>1){
      quimin=quimin-1
      Zopt=Z;Zopt[,-I2[ordre[1:quimin]]]=0
      if(plot){
         abline(v=quimin)
      }
   }else{
      Zopt=0*Z
   }
   if(!is.null(trunc)){
      trunc=min(trunc,p2)
      trunc=max(0,trunc)
      Zopt=Z;Zopt[,-I2[ordre[1:trunc]]]=0
      if(plot){
         abline(v=trunc,col="red")
      }
   }
   return(list(curve=curve,Zopt=Zopt))
}