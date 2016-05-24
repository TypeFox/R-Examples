krzanowski.test<-function(CA, CB, vecsA, vecsB, corr=FALSE, ...){

if(length(vecsA)!=length(vecsB)){stop("vecsA amd vecsB must be same length")}
if(length(vecsA)>c(dim(CA)[1]/2)){warning("component spaces must be ", floor(c(dim(CA)[1]/2)), " dimension","s"[floor(c(dim(CA)[1]/2))>1] ," or less")}

  if(corr){
    CA<-cov2cor(CA)
    CB<-cov2cor(CB)
  }

  CA<-as.matrix(eigen(CA)$vectors[,vecsA])
  CB<-as.matrix(eigen(CB)$vectors[,vecsB])

  N<-t(CA)%*%CB%*%t(CB)%*%CA

  eigval<-eigen(N)$values
  eigvec<-as.matrix(eigen(N)$vectors)
  bisector<-matrix(NA, dim(CA)[1], length(vecsA))
  I<-diag(dim(CA)[1])
  
  for(i in 1:length(vecsA)){ 
  bisector[,i]<-sqrt(2*(1+sqrt(eigval[i])))*((I+(1/sqrt(eigval[i]))*CB%*%t(CB))%*%CA%*%eigvec[,i])
  bisector[,i]<-bisector[,i]/sqrt(sum(bisector[,i]^2))
  }
  
  angles<-acos(eigval^0.5)/(pi/180)
  sumofS<-sum(eigval) 

  list(sumofS=sumofS, angles=angles, bisector=bisector)
}




 
