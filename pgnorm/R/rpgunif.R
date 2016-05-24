rpgunif <-
function(n,p){

# A function implemented by Steve Kalke

# Description: 
# Samples from the bivariate p-generalized uniform distribution

# Arguments: 
# p- a positiv constant (default: p=2)
# n- the number of random vectors to be simulated

if (missing(p)){p<-2}

Phi<-rpgangular(n,p)
U<-cbind(cos(Phi)/((abs(cos(Phi))^p+abs(sin(Phi))^p)^(1/p)),sin(Phi)/((abs(cos(Phi))^p+abs(sin(Phi))^p)^(1/p)))
return(U)
}
