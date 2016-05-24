data_simu <-
function(n,N,eta_star,q)
{
sigma_u=1
P=runif(N,0.1,0.5)
W=matrix(0,n,N)
for(j in 1:N)
{ W[,j]=rbinom(n,2,P[j])}
nb_comp_non_zero=q*N
sigma_e=sqrt(q*N*sigma_u^2*(1-eta_star)/eta_star)
b=sample(1:N,nb_comp_non_zero)
a1=sort(b)
u=rnorm(nb_comp_non_zero,0,sigma_u)
e=rnorm(n,0,sigma_e)
U=matrix(0,N) 
U[a1]=u
Z=scale(W,center=TRUE,scale=TRUE)
Y=Z%*%U + e
list(W=W,Y=Y)
}
