eltestwtinC<-function(x,wt,mu){
	pi=x;
	Lx=length(x);
	lam0=0;
 	if(any(is.na(x))) stop('NaNs');    ## if (sum(is.na(x))>0)stop('NaNs');	
 	re=.C('eltestwt',
		x=as.numeric(x),
		wt=as.numeric(wt),
		mu1=as.double(mu),   ##added as.double() 3/2015 by M Zhou
		Lx1=Lx,
		pi=as.numeric(pi),
		lamre=as.double(lam0)  ##added as.double() 3/2015 by M Zhou
	     )
	return(list(x=re$x, wt=re$wt, prob=re$pi, lam=re$lamre));
 }


el.test.wt <- function(x, wt, mu, usingC=TRUE) {
#x <- as.matrix(x)
#if( ncol(x) != 1 ) stop("x must be a vector") 

if( length(mu) != 1 ) stop("mu must be a scalar")

if(usingC){
return(eltestwtinC(x,wt,mu));
}else{
xmu <- x-mu
allw <- sum(wt)
BU <- 0.02*allw/max(abs(xmu))

lamfun <- function(lam,xmu,wt,allw) { sum(wt*xmu/(allw+lam*xmu)) }

if(lamfun(0,xmu,wt,allw) == 0) lam0 <- 0 
else {
 if( lamfun(0,xmu,wt,allw) > 0 ) {lo <- 0
                                up <- BU
                                while(lamfun(up,xmu,wt,allw)>0)
                                     up <- up + BU
                                 }
 else {up <- 0
      lo <- - BU
      while(lamfun(lo,xmu,wt,allw) < 0 )
           lo <- lo - BU
     }
 lam0 <- uniroot(lamfun,lower=lo,upper=up,tol=1e-9,xmu=xmu,wt=wt,allw=allw)$root
}
pi <- wt/(allw + lam0*xmu)
return(list(x=x, wt=wt, prob=pi, lam=lam0));
}
}
### add output of lam0, 5/2007 
