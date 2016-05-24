`hbpp.integrate1` <-
function(MU,V,A=1,RP=.1,...){
   1-pnorm(log10((1-RP)/RP)/A,mean=MU,sd=sqrt(V))
}

