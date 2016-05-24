tau2par.bvn=function(tau)
{ sin((pi*tau)/2) }

tau2par.cln=function(tau)
{ 2*tau/(1-tau) }

tau2par.cln180=function(tau)
{ 2*tau/(1-tau) }

tau2par.cln90=function(tau)
{ 2*tau/(1+tau) }

tau2par.cln270=function(tau)
{ 2*tau/(1+tau) }

tau2par.frk=function(tau){
  a<-1
  if(tau<0){
    a<- -1
    tau<- -tau}
  f = function(x) {
    x/(exp(x) - 1)}
  tauF=function(x) 1 - 4/x + 4/x^2 * integrate(f, lower = 0+.Machine$double.eps^0.5, upper = x)$value
  v<-uniroot(function(x) tau - tauF(x) ,lower=0,upper=500, tol = .Machine$double.eps^0.5)$root
  return(a*v)
}