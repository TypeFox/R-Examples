##A is the (already) projected & sphered data
##Projections must be sphered for these indices
ecdf.indices <- function(A, sphered = FALSE)
    {
        if(sphered == FALSE)
            {
                A <- princomp(A)$scores
            }
        
        s1 <- pnorm(A[,1])
        s2 <- pnorm(A[,2])
        
        n <- nrow(A)
        
        ecdf2 <- function(x)
            sum( (A[,1] <= x[1]) & (A[,2] <= x[2]) ) / n

        Fx <- sapply(1:n, FUN = function(i) ecdf2(A[i,]))
        Gx <- sapply(1:n, FUN = function(i) ecdf2(rev(A[i,])))
        
        K1 <- Fx - s1*s2
        K2 <- Fx - Gx

        return(c(CvM = sum(K1 ^ 2),
                 KS = max(abs(K1)),
                 D2 = sum(K2 ^ 2),
                 Dinf = max(abs(K2)) ))
    }
