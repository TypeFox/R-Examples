qpolicy_loss <-
function(q,mu,delta,lambda,theta,family,y.max=20,zt=TRUE){
	foo<-function(s){
	ppolicy_loss(s,mu,delta,lambda,theta,family,y.max,zt) -q
	}
    out<-uniroot(foo,lower=0,upper=(mu +4*mu*sqrt(delta))*2*lambda)$root
    return(out)
}
