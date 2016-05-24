aip <-
function(p_prior, sd=1.34){
# Singletons calculations in order that E(psi) = p_prior
fu<- function(s,p_prior){sapply(s,FUN= function(s){integrate(function(a){psip(s,a)*fp(a, sd)},-Inf,Inf)$value - p_prior})}
sgl<- vector()
for (i in 1:length(p_prior)){
fu_i<- function(s){fu(s,p_prior[i])}
sgl[i]<- uniroot(fu_i,c(0,1))$root
}
sgl
}
