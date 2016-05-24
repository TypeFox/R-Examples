###krylow.pls<-function(D,d,dim)
###{
###R=d;  Sxxsxy=R;
###if (dim>=2)
###for (i in 2:dim)
###{
###Sxxsxy=D %*% Sxxsxy  ;
###R=cbind(R,Sxxsxy);
###}
###beta= R %*% solve(t(R) %*% D %*% R) %*% t(R) %*% d
###beta=  beta;
###return(list(beta=beta))
###}

### Anders Gorst Rasmussen's code
krylow.pls <- function(D,d,dim = 1) {
	r <- d
	p <- r
	rsold <- drop(t(r) %*% r) 
	x <- rep(0,length(d)) 
	for(k in 1:dim){
		Ap <- D %*% p;
		alpha <-drop(rsold / (t(p) %*% Ap)); x <- x + alpha * p;
		r <- r - alpha * Ap;
		rsnew <- drop(t(r) %*% r);
		p <- r + rsnew / rsold * p;
		rsold <-rsnew;
	}
	return(x) 
}
