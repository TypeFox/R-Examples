ail2 <-
function(p_prior){
# Singletons calculations in order that E(psi) = p_prior
fu<- function(s,p_prior){sapply(s,FUN= function(s){integrate(function(a){psil(s,a)*fl(a)},0,100)$value - p_prior})}
sgl<- vector()
for (i in 1:length(p_prior)){
fu_i<- function(s){fu(s,p_prior[i])}
sgl[i]<- uniroot(fu_i,c(-105,5))$root
}
sgl
}
