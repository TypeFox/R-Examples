ardec.periodic <-
function(x,per,tol=0.95){

# if(frequency(x)!=12){stop("monthly time series required")}
# updated 29 Apr 2013

fit=ardec.lm(x)

comp=ardec(x,fit$coefficients)

if(any(comp$period > (per-tol) & comp$period < (per+tol))) {
                candidates=which(comp$period > (per-tol) & comp$period < (per+tol))
                lper=candidates[which.max(comp$modulus[candidates])]
                l=comp$period[lper]
                m=comp$modulus[lper]
                gt=Re(comp$comps[lper,]+comp$comps[lper+1,])   }

return(list(period=l,modulus=m,component=gt))  

}
