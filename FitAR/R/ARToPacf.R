`ARToPacf` <-
function(phi){
phik=phi
L=length(phi)
if(L==0) return(0)
pi=numeric(L)
for (k in 1:L){
    LL=L+1-k
    a <- phik[LL]
    pi[L+1-k] <- a
    phikp1 <- phik[-LL]
    if(is.na(a) || abs(a)==1)
        stop("transformation is not defined, partial correlation = 1")
    phik <- (phikp1+a*rev(phikp1))/(1-a^2)
    }
pi
}
