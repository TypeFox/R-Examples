# Input

# y input data
# M maximum template length
# r matching tolerance

# Output

# se standard error estimates for m=0,1,...,M-1
# e  sample entropy estimates for m=0,1,...,M-1

.packageName <- 'mousetrack'

sampense <- function(y, M, r){

    res = makerun(y,M,r)
    
    F1 = res$F1; F2 = res$F2; R1 = res$R1; R2 = res$R2
    
    F = F1 + F2
    n = length(y)
    dd = 1
    K0 = colSums(F*(F-1))
    K = matrix(0, nrow = M, ncol = M) 
    
    for (m in 1:M){
        kbuild = vector()
        
        for (d in 1:min(m+1,M)){
            
            i1 = (d+1):n
            i2 = i1-d     
            nm1 = F1[i1,m]
            nm2 = F2[i1,m]
            nm3 = F1[i2,m]
            nm4 = F2[i2,m]
            
            nm1 = nm1 - vector("numeric",length(i1) ) # sum(R2(i2,1:(dd-1)) >= m )
            
            fill2 = vector("numeric", length = length(i1) )
            ind2 = which(R2[i1, 1:(2*d)] >= m, arr.ind = TRUE )
            
            for (f2 in unique(ind2[,1]) ){
                ind = which(ind2[, 1] == f2)
                fill2[f2] = length(ind) 
            }
            
            nm2 = nm2 - fill2 
            
            if (d > 1){
                
                fill3 = vector("numeric", length = length(i2) )
                ind3 = which(R1[i2, 1:(2*d-1)] >= m, arr.ind = TRUE )
                
                for (f3 in unique(ind3[,1]) ){
                    ind = which(ind3[, 1] == f3)
                    fill3[f3] = length(ind) 
                }
                
                nm3 = nm3 - fill3 

            } else {
                
                fill3 = vector("numeric", length = length(i2) )
                ind3 = which(R1[i2, 1:(2*d-1)] >= m)
                fill3[ind3] = 1   
                nm3 = nm3 - fill3
                
            }
            
            nm4 = nm4 - vector("numeric",length(i2) ) # sum(R2(i2,1:(dd-1)) >= m )  
            K[d,m] = 2*sum((nm1+nm2)*(nm3+nm4))
        }
    }

    K = rbind(K0, K, deparse.level = 0)

    n1 = n2 =  vector("numeric", length = M)

    n1[1] = n*(n-1)*(n-2)
    
    for (m in 1:(M-1)){
        n1[m+1] = sum(K[1:(m+1),m])
    }
    
    for (m in 1:M){
        n2[m] = sum(K[1:m,m])
    }
    
    A = colSums(F)/2
    N = n*(n-1)/2
    B = c(N, A[1:(M-1)])
    p = A/B;
    e = -log(p);
    
    vp = p*(1-p)/B + (n2-n1*p^2)/(B^2) ## matlab function has a weird max
    sp = sqrt(vp)
    se = sp/p
    
    return(list (se = se, e = e) )

}

    
