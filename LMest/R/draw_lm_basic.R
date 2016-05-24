draw_lm_basic <-
function(piv,Pi,Psi,n){

#        [lk,piv,Pi,Psi,np,aic,bic,lkv] = draw_lm_basic(piv,Pi,Psi,n)
#
# Draw a sample of size n from a Basic Latent Markov model with parameter piv, Pi and Psi

# Preliminaries
	k = length(piv)
	dd = dim(Psi)
	c = dim(Psi)[1]
	TT = dim(Pi)[3]
	if(length(dd)>2) r = dd[3]
	else r = 1
# For each subject
    Y = matrix(0,n,TT*r)
    cat("------------|\n")
    cat(" sample unit|\n")
    cat("------------|\n")
    for(i in 1:n){
    	if(i/1000==floor(i/1000)) cat(sprintf("%11g",i),"\n",sep=" | ")
    	u = k+1-sum(runif(1)<cumsum(piv))
    	ind = 0
    	for(j in 1:r){
    		ind = ind+1
	    	Y[i,ind] = c-sum(runif(1)<cumsum(Psi[,u,j]))		
    	}
   		for(t in 2:TT){
    		u = k+1-sum(runif(1)<cumsum(Pi[u,,t]))
	    	for(j in 1:r){
    			ind = ind+1
    			Y[i,ind] = c-sum(runif(1)<cumsum(Psi[,u,j]))
    		}	
   		}
    }
    cat("------------|\n")
	out = aggr_data(Y)
    S = out$data_dis; yv = out$freq
    S = array(t(S),c(r,TT,length(yv)))
    S = aperm(S)
    if(r==1) S = S[,,1] 
    out = list(Y=Y,S=S,yv=yv)
}
