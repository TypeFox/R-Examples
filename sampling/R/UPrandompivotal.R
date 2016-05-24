"UPrandompivotal" <-
function(pik,eps=1e-6)
{
if(any(is.na(pik))) stop("there are missing values in the pik vector")
N=length(pik)
v=sample(c(N),N)
s=numeric(N)
s[v]=UPpivotal(pik[v],eps)
s
}

