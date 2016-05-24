MLM.beadarray <-
function(sig,stderr,nbeads,groups,var.equal = FALSE,max.iteration=20,epsilon=1e-6,method="REML"){
 if(method!="REML") {method="ML"}
 K <- length(groups)
 n <- sapply(groups,length)
 theta<-tau2<-sigma2<-NULL; previous.theta <- matrix(rep(0,K*nrow(sig)),ncol=K)
 for(k in 1:K){
   g<-groups[[k]]
   theta  <-cbind( theta,rowMeans(sig[,g])) 
   tau2   <-c(  tau2,list(pmax((rowSums(sig[,g]^2) - rowSums(sig[,g])^2/n[k]) / (n[k]-1)- rowSums(stderr[,g]^2)/n[k],0))) 
   sigma2 <-c(sigma2,list(stderr[,g]^2*nbeads[,g])) 
 }
 if(var.equal == TRUE) {tau2common <-  rowMeans(sapply(tau2,function(x){x}))}  
 for(j in 1:max.iteration){
   w <- A <- B <- C <- D <- E<-NULL 
     for(k in 1:K){
      g<-groups[[k]]
      if(var.equal == FALSE) { 
         w <-c(w,list(1/(sigma2[[k]]/nbeads[,g]+tau2[[k]]))) 
      }
      if(var.equal == TRUE) { 
         w <-c(w,list(1/(sigma2[[k]]/nbeads[,g]+tau2common)))
      }
      B <-c(B,list((w[[k]]/nbeads[,g])^2/((w[[k]]/nbeads[,g])^2 + (nbeads[,g]-1)/(stderr[,g]^4*nbeads[,g]^2)))) 
      A <-cbind(A,rowSums(w[[k]]^2)-rowSums(B[[k]]*w[[k]]^2/nbeads[,g]))
      C <-cbind(C,rowSums( (1-B[[k]]/nbeads[,g])*(w[[k]]^2*(sig[,g]-theta[,k])^2-w[[k]])+B[[k]]*((nbeads[,g]-1)*(sigma2[[k]]-stderr[,g]^2*nbeads[,g])/sigma2[[k]]^2   )))
      D <-c(D,list((nbeads[,g]/w[[k]])^2*((w[[k]]^2*(sig[,g]-theta[,k])^2-w[[k]])/nbeads[,g]-(nbeads[,g]-1)*(sigma2[[k]]-stderr[,g]^2*nbeads[,g])/sigma2[[k]]^2)))
      if(method=="REML") {E <-cbind(E,rowSums(w[[k]]^2)/rowSums(w[[k]]))}
     } 
   if(var.equal == FALSE) { 
     for(k in 1:K){
      if(method=="REML") {tau2[[k]] <- pmax(tau2[[k]] + C[,k]/A[,k] + tau2[[k]]^2*E[,k]/n[k],0)} else {tau2[[k]] <- pmax(tau2[[k]] + C[,k]/A[,k],0)}
      sigma2[[k]] <- sigma2[[k]] +B[[k]]*D[[k]]- B[[k]]*C[,k]/A[,k]
     }  
   }
   if(var.equal == TRUE) { 
    totalA <- rowSums(A)
    totalC <- rowSums(C)
    if(method=="REML") {totalE <- rowSums(E)}
    if(method=="REML") {tau2common <- pmax(tau2common + totalC/totalA + tau2common^2*totalE/sum(n),0)} else {tau2common <- pmax(tau2common + totalC/totalA,0)}
    sigma2[[k]] <- sigma2[[k]] +B[[k]]*D[[k]]- B[[k]]*totalC/totalA
   }
   w<-theta<- NULL
   for(k in 1:K){
     g<-groups[[k]]
     sigma2[[k]][sigma2[[k]]<0]<-0
     if(var.equal == FALSE) { w <-c(w,list(1/(sigma2[[k]]/nbeads[,g]+tau2[[k]]))) }
     if(var.equal == TRUE) { w <-c(w,list(1/(sigma2[[k]]/nbeads[,g]+tau2common))) }
     theta <-cbind(theta,rowSums(w[[k]]*sig[,g])/rowSums(w[[k]])) 
   }
  convergence.check <- rowSums(abs(theta-previous.theta) > epsilon)>0
  print(paste("Number of converged transcripts: ",sum(!convergence.check),". Now ",sum(convergence.check)," remaining...",sep=""))
  if(j>1 & sum(convergence.check) == 0 ) {break}
  previous.theta <- theta
  } # for j
  v<-NULL 
  for(k in 1:K){v <- cbind(v,1/rowSums(w[[k]]))}
  if(var.equal == FALSE) { res<-(c(theta=list(theta),tau2=tau2,sigma2=sigma2))}
  if(var.equal == TRUE) { res<-(c(theta=list(theta),tau2=list(tau2common),sigma2=sigma2))}
  ###print("Computing the Wald test statistics...")
  if(K == 2) {
    t.statistics=(theta[,2]-theta[,1])/sqrt(v[,1]+v[,2])
    res<- c(res,t.statistics=list(t.statistics),method=method,convergence=list(!convergence.check))
  } 
  if(K > 2) {
   f.statistics <- 0
   for(k in 1:(K-1)){  for(h in (k+1):K){   f.statistics <-   f.statistics + (theta[,k]-theta[,h])^2/(v[,k]*v[,h]) } }
   const <- rowSums(1/v)
   res<- c(res,f.statistics=list((1/(K-1))*f.statistics/const),method=method,convergence=list(!convergence.check))
  } 
  return(as.data.frame(res))
}
