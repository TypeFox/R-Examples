ChangeX <- function(N, geno, covariates, num.test)
{
    # geno is matrix or data.frame
    # covariates is NULL or matrix

    if (is.null(covariates))
    {
        Len <- (1 + length(covariates) + ncol(geno)) * num.test
    }else
    {
        Len <- (1 + ncol(covariates) + ncol(geno)) * num.test
    }

    x.mat <- matrix(data=0, nrow=num.test*N, ncol=Len)

    change.id <- matrix(data=1:(num.test*N), nrow=N, ncol=num.test)

    dex.rec <- change.id[,1]
    dex.add <- change.id[,2]
    dex.dom <- change.id[,3]
    change.id <- as.vector(t(change.id))

    geno.rec <- geno
    geno.add <- geno
    geno.dom <- geno
    geno.rec[geno.rec==1] <- 0
    geno.rec[geno.rec==2] <- 1
    geno.dom[geno.dom==2] <- 1

    mat.rec <- cbind(1, covariates, geno.rec)
    mat.add <- cbind(1, covariates, geno.add)
    mat.dom <- cbind(1, covariates, geno.dom)
    
    k <- ncol(mat.rec)
    list.1 <- 1:k
     
    x.mat[dex.rec, list.1] <- mat.rec
    x.mat[dex.add, list.1+k] <- mat.add
    x.mat[dex.dom, list.1+k*2] <- mat.dom

    x.mat[change.id,]
}
