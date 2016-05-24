`Kmulti.ls` <-
function (X, I, J, r = NULL, corre = "isotropic") 
{
    n1 = sum(I)
    n2 = sum(J)
    cosa12 <- Kmulti(X, I, J, r, correction = corre)
    cosa21 <- Kmulti(X, J, I, r, correction = corre)
    K12ls = ((n2 * cosa12[[3]]) + (n1 * cosa21[[3]]))/(n1 + n2)
    cosa12[[3]] = K12ls
    return(cosa12)
}

