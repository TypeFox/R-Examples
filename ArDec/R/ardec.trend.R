ardec.trend <-
function(x){

options(warn=-1)

fit=ardec.lm(x)


comp=ardec(x,fit$coefficients)

if(any(comp$period==Inf)){warning("no trend component")}


if(any(comp$period ==Inf)){
                l=comp$period[which(match(comp$period,Inf)==1)[1]]
                m=comp$modulus[which(match(comp$period,Inf)==1)[1]]
                gt=Re(comp$comps[which(match(comp$period,Inf)==1 )[1],])

 }


return(list(modulus=m,trend=gt))  





}
