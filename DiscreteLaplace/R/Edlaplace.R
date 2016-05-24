Edlaplace <-
function(p,q)
{
E1<-1/(1-p)-1/(1-q)
E1a<-(q*(1-p)^2+p*(1-q)^2)/((1-q*p)*(1-q)*(1-p))
V<-1/(1-p)^2/(1-q)^2*((q*(1-p)^3*(1+q)+p*(1-q)^3*(1+p))/(1-p*q)-(p-q)^2)
list(E1=E1,E1a=E1a,V=V)
}

