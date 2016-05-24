## translated to R by Moreno I. Coco (moreno.cocoi@gmail.com)
## GPL < 2

## Input

## y input data
## M maximum template length
## r matching tolerance

## Output

## e sample entropy estimates for m=0,1,...,M-1
## A number of matches for m=1,...,M
## B number of matches for m=0,...,M-1 excluding last point

.packageName <- 'mousetrack'

sampenc <- function(y, M, r){

    n = length(y)

    run = lastrun = vector(mode = "numeric", length = n)
    A = B =  vector(mode = "numeric", length = M)

    for (i in 1:(n-1)){

        nj = n-i
        y1 = y[i]
    
        for (jj in 1:nj){
            j = jj+i
        
            if (abs(y[j] - y1) < r){
                
                run[jj] = lastrun[jj] + 1
                M1 = min(M, run[jj])

                for (m in 1:M1){           
                    A[m] = A[m] + 1
                
                    if (j < n){
                        B[m] = B[m] + 1
                    }
                }
            } else {
                run[jj] = 0
            }
        }
        
        for (j in 1:nj){
            lastrun[j] = run[j]
        }
    }


N = n*(n-1)/2;
B = c(N, B[1:(M-1)])
p = A/B;
e = -log(p);

return ( cbind(e, A, B) )

}
