"weightslogl.calc.w" <- 
function(p,fjk,weights){
# amended to rescale posterior probabilities
# R.E.D. 15th Feb 2006
# p is a vector of length K containing the mixture proportions
# fjk is a JXK matrix of log densities
logpf <- t(apply(fjk,1,"+",log(p)))
Mi <- apply(logpf,1,max)
Sik <- logpf-Mi
ifelse(Sik < -760, 0, Sik) 
eSik <- exp(Sik)
SeSik <- as.vector(apply(eSik,1,sum))
w <- eSik/SeSik
ML.dev <- -2*sum(weights*log(apply(eSik,1,sum)))-2*sum(weights*Mi)
list(w=w, ML.dev=ML.dev)
} 



