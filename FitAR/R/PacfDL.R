`PacfDL` <-
function(c, LinearPredictor=FALSE){
L<-length(c)-1
d<-c[-1]
if (L==0) {
    phik<-numeric(0)
    vk<-c[1]
    pi<-numeric(0)
}
else {
    phik<-c[2]/c[1]
    pi<-numeric(L)
    pi[1]<-phik
    vk <- c[1]*(1 - pi[1]^2)
}
if (L>1){
   for (k in 2:L) {
        vkm1 <- vk
        phikm1 <- phik
        a <- sum(c(1,-phikm1)*rev(d[1:k]))/vk
        phik <- c(phikm1-a*rev(phikm1),a)
        vk <- vkm1*(1-a^2)
        pi[k] <- a
    }
}
if (!LinearPredictor)
    list(Pacf=pi)
else
    list(Pacf=pi, ARCoefficients=phik, ResidualVariance=vk)
}

