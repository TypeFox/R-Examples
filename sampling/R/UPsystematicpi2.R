UPsystematicpi2<-function(pik) 
{
    n=sum(pik)
    n=.as_int(n)
    pik1 = pik[pik > 0 & pik < 1]
    N = length(pik1)
    Vk = cumsum(pik1)
    Vk1=Vk%%1 
    if(Vk1[N]!=0) Vk1[N]=0
    r = c(sort(Vk1), 1)
    cent = (r[1:N] + r[2:(N + 1)])/2
    p = r[2:(N + 1)] - r[1:N]
    A = matrix(c(0, Vk), nrow = N + 1, ncol = N) - 
    t(matrix(cent,nrow = N, ncol = N + 1))
    A = A%%1
    M = matrix(as.integer(A[1:N, ] > A[2:(N + 1), ]), N, N)
    pi21 = M %*% diag(p) %*% t(M)
    pi2 = pik %*% t(pik)
    pi2[pik > 0 & pik < 1, pik > 0 & pik < 1] = pi21
    pi2
}


