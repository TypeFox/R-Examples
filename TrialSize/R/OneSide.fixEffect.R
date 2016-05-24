OneSide.fixEffect <-
function(m,m1,delta,a1,r1,fdr){

a2=1-a1
alpha_star=r1*fdr/((m-m1)*(1-fdr))
beta_star=1-r1/m1
n=floor((qnorm(1-alpha_star)+qnorm(1-beta_star))^2/(a1*a2*delta^2))+1
return(n)
}
