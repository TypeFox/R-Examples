mNNinfo2b <- 
function (n, R, Q) 
{
    N <- sum(n)
    SP <- n # cambiamos n por SP (numero individuos por especie) para que no haya dos variables con el mismo nombre (N y n) en Fortran
    k <- length(n)
    l <- names(n)
    EN <- matrix(0, nrow = k, ncol = k)
    VN <- VarN <- matrix(0, nrow = k * k, ncol = k * k)
    for (i in 1:k) {
        for (j in 1:k) {
            EN[i, j] <- n[i] * (n[j] - (i == j))/(N - 1)
        }
    }
   
                    
    ans <- .Fortran('dixonloop', k=as.integer(k), N=as.double(N),
                    R=as.double(R), Q= as.double(Q), SP=as.double(SP),
                    VN=as.integer(VN), VarN=as.double(VarN), EN=as.double(EN),
                    PACKAGE="dixon")                    

    VN <- matrix(ans$VN,nrow = k * k, ncol = k * k)
    VarN <- matrix(ans$VarN,nrow = k * k, ncol = k * k)
    EN <- matrix(ans$EN,nrow = k, ncol = k)
    
    v <- as.vector(t(outer(l, l, paste, sep = ""))) # ojo "sep". igual con sep="." mejor
    dimnames(VN) <- dimnames(VarN) <- list(v, v)
    list(EN = EN, VN = VN, VarN = VarN)
}
