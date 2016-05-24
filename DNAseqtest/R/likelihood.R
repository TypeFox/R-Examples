likelihood <-
function(thetast,fobs,merge2){
fobss=fobs;thetass=thetast;merge22=merge2
fst<-gn(thetass,merge22)
if(min(fst)<=0){
return(lik<-4)
}
lik<-sum(fobss*log(fobss/fst))
lik
}
