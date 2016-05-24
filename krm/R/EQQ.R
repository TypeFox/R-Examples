#Compute the second non-central moment of a score.
E.QQ = function(K, mu, do.C=TRUE) {
    
    if (!is.matrix(K)) K=as.matrix(K)
    if (nrow(K)!=ncol(K)) stop("K is not a square matrix.")
    n=nrow(K)
    if (length(mu)!=n) stop("mu is not of length n.")
    
    mom=0
    
    if (do.C) {
    
        if (!is.double(K)) K <- as.double(K)    
        aux=.C("score_var", "_K" = K, "_n"=as.integer(n), "_mu"=as.double(mu), "_mom"=as.double(mom))
        mom=aux$"_mom"
    
    } else {
    
        mom = mom + sum((mu*(1-mu)^4 + (1-mu)*mu^4) * diag(K)^2)

        for (i in 1:n) {
            for (k in 1:n) {                        
                if (k==i) next                
                mom = mom + mu[i]*(1-mu[i]) * mu[k]*(1-mu[k]) * K[i,i] * K[k,k]
        }}

        for (i in 1:n) {
            for (j in 1:n) {                        
                if (j==i) next                
                mom = mom + mu[i]*(1-mu[i]) * mu[j]*(1-mu[j]) * K[i,j] * K[i,j] * 2
        }}

    }
    
    mom 
    
}
## testing
#n=1e1
#K=diag(n)
## make it a block diagonal
#for (i in 1:n) {
#    if (i \%\% 2 ==1) K[i,i+1]=1 else K[i,i-1]=1
#}
#
#E.QQ (K, mu=rep(1/2,n), do.C=TRUE)
#E.QQ (K, mu=rep(1/2,n), do.C=FALSE)
