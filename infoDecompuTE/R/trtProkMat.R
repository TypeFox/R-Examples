	#pre- and post-multiply NTginvATN by block projection matrices produces the effFactors for each stratum  
trtProkMat <-            
function(z, T, N, Rep){
if(!is.matrix(z)) return(z)

nEffect = length(T)
PNTginvATNP = effFactors = vector("list", nEffect)

names(PNTginvATNP) = names(effFactors) = names(T)

PNTginvATNP[[1]] = z %*% N %*% T[[1]] %*% invInfMat(C=z, N=N, T=T[[1]]) %*% T[[1]] %*% t(N) %*% t(z) 
    
if(!all(PNTginvATNP[[1]] <1e-6)){    
  effFactors[[1]] = vector("list", nEffect)
  names(effFactors[[1]]) = names(T)
  for(i in 1:nEffect){
     r.adjust = ginv(sqrt(diag(Rep[,i])))
    #eigenvalues of the information matrix
    va = Re(eigen(r.adjust %*% T[[i]] %*% t(N) %*% PNTginvATNP[[1]] %*% N %*% T[[i]] %*% r.adjust)$va)
    #va = Re(eigen(T[[i]] %*% t(N) %*% PNTginvATNP[[1]] %*% N %*% T[[i]])$va)
   
    #harmonic means of the canonical efficiency factors to give the average efficiency factor
    effFactors[[1]][[i]] = 1/mean(1/va[which(va>1e-6)]) 
  } 
}

newZ = (z %*% t(z)) - PNTginvATNP[[1]]

if(nEffect !=1){            
  for(i in 2:nEffect){                      
    PNTginvATNP[[i]] = newZ %*% N %*% t(T[[i]]) %*% invInfMat(C=newZ, N=N, T=T[[i]]) %*% T[[i]] %*% t(N) %*% t(newZ)
            
    if(all(PNTginvATNP[[i]] <1e-6)) next
    
    effFactors[[i]] = vector("list", nEffect)
    names(effFactors[[i]]) = names(T)
    for(j in 1:nEffect){
       r.adjust = ginv(sqrt(diag(Rep[,j])))
      va = Re(eigen(r.adjust %*% T[[j]] %*% t(N) %*% PNTginvATNP[[i]] %*% N %*% T[[j]]%*% r.adjust)$va)
      #va = Re(eigen(T[[j]] %*% t(N) %*% PNTginvATNP[[i]] %*% N %*% T[[j]])$va)
     
      #effFactors[[i]][[j]] = 1/mean(Rep[names(T[j])]/va[which(va>1e-6)])
      effFactors[[i]][[j]] = 1/mean(1/va[which(va>1e-6)])
    
    }
    newZ = (newZ %*% t(newZ)) - PNTginvATNP[[i]]  
   
    }
} 

return(effFactors) 
}
