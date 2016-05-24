calcR2 <- function(matrix_pol,populations,thetaT,S){

########################
# when include_unknown=0
# Haploid
# Ohne Outgroup, da in print_output eh nich dabei !
########################

npops <- length(populations)
R2    <- rep(NaN,npops)
nn    <- paste("pop",1:npops)
 
for(x in 1:length(populations)){
   mat_rix <- matrix_pol[populations[[x]],]
   
 if(is.matrix(mat_rix)){
     unic <- numeric(dim(mat_rix)[1])
     freq <- unic 
     erg <- apply(mat_rix,2,ap_calcR2,unic,freq)
 }else{
     unic <- numeric(length(mat_rix))
     freq <- unic 
     erg <- sapply(mat_rix,ap_calcR2,unic,freq)
 }
#if(is.matrix(mat_rix)){erg <- apply(mat_rix,2,ap_calcR2)}
#if(is.vector(mat_rix)){erg <- sapply(mat_rix,ap_calcR2)}

 ifelse(is.matrix(erg), p_unic <- rowSums(erg),p_unic <- sum(erg))
 
#if(is.matrix(erg)){ p_unic <- rowSums(erg)}
#if(is.vector(erg)){ p_unic <- sum(erg)}

R2[x] <- R2(p_unic,thetaT[x],length(populations[[x]]),S[x])
}
names(R2) <- nn

return(R2)

}

#APPLY FUNCTION
ap_calcR2 <- function(vek,unic,freq){

   freq[2] <-sum(vek==0,na.rm=TRUE)  
   freq[3] <-sum(vek==1,na.rm=TRUE)  
   freq[4] <-sum(is.na(vek)) 
   freq[1] <- freq[2] +freq[3] 
   
# UNIC NUMBER OF SINGELTONS IN EACH SEQUENCE
if(freq[1]){
  if(freq[2]==1){
     id <- vek==0
     unic[id] <- 1
  }
  if(freq[3]==1){
    id <-  vek==1
    unic[id] <- 1
  } 

}#END if(freq[1])


return(unic)
}
##### Function R2

R2 <- function(unic,pi,sample_size,S)
{
    # pi is thetaT
    # sample_size -> Groesse der Population 

    sm2 <- 0
    
    if(S == 0 || sample_size == 0){ return(NaN)}
    
    sm2 <- sum((unic-pi/2)^2)
    
    #for (i in 1:sample_size){
		#sm2 <- sm2 + (unic[i] - pi/2)*(unic[i] - pi/2)
    #}
    
    sm2 <- sqrt(sm2/(sample_size))/S

    if(is.na(sm2)){return(NA)}
	
    if (abs(sm2) < 1.0e-15){
		sm2 <- 0
    }	
 return (sm2)
}
