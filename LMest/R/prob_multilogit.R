prob_multilogit <-
function(Xdis,be,label,fort=TRUE,der=FALSE){

# der = T if the derivative is required
	 k = dim(Xdis)[1]
     ndis = max(label)
     n = length(label)
     ncov = length(be)
     Pdis = matrix(0,ndis,k); P = matrix(0,n,k)
     if(fort==F){
         for(i in 1:ndis){
            pdis = exp(Xdis[,,i]%*%be); pdis = pdis/sum(pdis)
            Pdis[i,] = pdis
     		mul = sum(label==i)
     		P[label==i,] = rep(1,mul)%o%pdis
          }
     }else{
	     out = .Fortran("prob_multilogit",Xdis=Xdis,be=be,label=as.integer(label),Pdis=Pdis,P=P,k=as.integer(k),ndis=as.integer(ndis),
	     ns=as.integer(n),ncov=as.integer(ncov))
	     Pdis = out$Pdis; P = out$P
     }
     if(der){
        lbe = length(be)	
        dPdis = array(0,c(ndis,k,lbe)); dP = array(0,c(n,k,lbe))
		for(i in 1:ndis){
            pdis = Pdis[i,]
            Op = diag(pdis)-pdis%o%pdis
            dPdis[i,,] = Op%*%Xdis[,,i]
            ind = which(label==i)
     		for(j in ind) dP[j,,] = dPdis[i,,]
          }
     }else{
     	dPdis = NULL; dP = NULL
     }
     out = list(Pdis=Pdis,P=P,dPdis=dPdis,dP=dP)
     return(out)
}
