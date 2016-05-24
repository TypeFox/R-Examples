############################################################
# Input

# y input data
# M maximum template length
# r matching tolerance

# Output

# F1 matches with future points
# R1 runs with future points
# F2 matches with past points
# R2 runs with past points

.packageName <- 'mousetrack'

makerun <- function(y, M, r){

    n = length(y)
    run1 = vector("numeric", length = n)
    MM = 2*M

    R1 = R2 = matrix(0, nrow = n, ncol = MM)
    F = F1 = matrix(0, nrow = n, ncol = M)
    
    for (i in 1:(n-1)){
        j = (i+1):n
        
        match = abs(y[j]-y[i]) < r 
        k = which(match == TRUE)
        nj = length(j)
        run = vector("numeric", length = length(j))
        run[k] = run1[k]+1
        
        for (m in 1:M){
            k = which(run >= m)      
            nm = length(k)
            F1[i,m] = nm
            F[i,m] = F[i,m] + nm
            F[i+k,m] = F[i+k,m] + 1 
        }
        
        nj = min(nj, MM)
        k = (1:nj)
        R1[i,k] = run[k] ## maybe here issues on how to fill this in
        run1 = run   
    }
    
    for (i in 2:n){
        nj = min(MM, i-1)
        
        for (j in 1:nj){
        R2[i,j] = R1[i-j, j]
    }
    }
    
    F2 = F-F1

    return(list (F1 = F1, F2 = F2, R1 = R1, R2 = R2) )

}
