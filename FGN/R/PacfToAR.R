`PacfToAR` <-
function(zeta){
L=length(zeta)
if (L==0) return(numeric(0))
if (L==1) return(zeta)
phik=zeta[1]
for (k in 2:L){
    phikm1=phik
    phik=c(phikm1-zeta[k]*rev(phikm1),zeta[k])
    }
phik
}

