Mstep<-function(Z=Z,X=X,sigma_IR=sigma_IR,Ir=Ir){
   alpha=hatB(Z = Z,X =X )
   for (j in Ir){
      sigma_IR[j]=sd(X[,j]-X%*%alpha[-1,j])
   }
   return(list(alpha=alpha,sigma_IR=sigma_IR))#en C on fera un void
}