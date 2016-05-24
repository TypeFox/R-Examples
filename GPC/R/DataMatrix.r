DataMatrix <- function(A,Des,InputDistrib,pmaxi,PCSpace){
  #print("Entering DataMatrix")
  DesLength=nrow(Des)
  Pval=c()
  if(PCSpace=="Physic"){
    for(i in 1 : ncol(A)){
      if(InputDistrib[i]=="Gaussian"){
        PHerm = hermite.he.polynomials(pmaxi, normalized=FALSE)
        facto=sqrt(matrix(t(matrix(rep(factorial(A[,i]),DesLength),nrow(A),DesLength)),1,nrow(A)*DesLength))
        Pval = rbind(Pval,unlist(polynomial.values(PHerm[(A+1)[,i]],Des[,i]))/facto)
      }
      else if(InputDistrib[i]=="Uniform") {
        PLeg = legendre.polynomials(pmaxi, normalized=FALSE)
        facto=sqrt(matrix(t(matrix(rep(2*(A[,i])+1,DesLength),nrow(A),DesLength)),1,nrow(A)*DesLength))
        Pval = rbind(Pval,unlist(polynomial.values(PLeg[(A+1)[,i]]*Des[,i]))*facto)
      }
      else if(InputDistrib[i]=="Gamma") {
        PLag = laguerre.polynomials(pmaxi, normalized=FALSE)
        facto=sqrt(matrix(t(matrix(rep(2*(A[,i])+1,DesLength),nrow(A),DesLength)),1,nrow(A)*DesLength))#Mettre ? jour
        Pval = rbind(Pval,unlist(polynomial.values(PLag[(A+1)[,i]],2*Des[,i]-1))*facto)
      }
      else if(InputDistrib[i]=="Beta") {
        PJac = jacobi.g.polynomials(pmaxi, normalized=FALSE)
        facto=sqrt(matrix(t(matrix(rep(2*(A[,i])+1,DesLength),nrow(A),DesLength)),1,nrow(A)*DesLength))#Mettre ? jour
        Pval = rbind(Pval,unlist(polynomial.values(PJac[(A+1)[,i]],2*Des[,i]-1))*facto)
      }
    }
  }
  else if(PCSpace=="Uniform"){
    for(i in 1 : ncol(A)){
      PLeg = legendre.polynomials(pmaxi, normalized=FALSE)
      facto=sqrt(matrix(t(matrix(rep(2*(A[,i])+1,DesLength),nrow(A),DesLength)),1,nrow(A)*DesLength))
      Pval = rbind(Pval,unlist(polynomial.values(PLeg[(A+1)[,i]],Des[,i]))*facto)
    }
  }
  if(PCSpace=="Gaussian"){ 
    for(i in 1 : ncol(A)){
      PHerm = hermite.he.polynomials(pmaxi, normalized=FALSE)
      facto=sqrt(matrix(t(matrix(rep(factorial(A[,i]),DesLength),nrow(A),DesLength)),1,nrow(A)*DesLength))
      Pval = rbind(Pval,unlist(polynomial.values(PHerm[(A+1)[,i]],Des[,i]))/facto)
    }
  }
  NPval = apply(t(Pval),1,FUN="prod")
  #print("Calculating poly matrix")
  M = matrix(NPval,DesLength,nrow(A))
  return(M)  
}