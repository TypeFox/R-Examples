bd.groups.optim<-function(phy,S,xcut=0,lambda=0,mu=0,survival=1) 
{
    if (class(phy) != "phylo") 
        stop("object \"phy\" is not of class \"phylo\"")
    if (!is.null(names(S))) {
        if (all(names(S) %in% phy$tip.label)) 
            S <- S[phy$tip.label]
        else warning("the names of argument \"S\" and the names of the tip labels\ndid not match: the former were ignored in the analysis.")
    }
    N <- length(S)
    x <- sort(branching.times(phy),decreasing=TRUE)
    trm.br <- phy$edge.length[phy$edge[, 2] <= N]
    dev <- function(par) {
    	a<-par[1]
    	r<-par[2]
        l <- r/(1-a)
 		m <- l*a  #r=l-m, a=m/l
 		if (l<=0 || m>l || m<0) {res <- 10^10} else {
 		p0m <- function(t,l,m){ (1-exp(-(l-m)*t))/(l-m*exp(-(l-m)*t)) }  #p0(t)/mu
 		p1 <- function(t,l,m){ (l-m)^2 * exp(-(l-m)*t)/(l-m*exp(-(l-m)*t))^2 }
 		if (xcut>0){trm.br<-trm.br*0+xcut}
 		#print(trm.br)
 		res<- -( (N-2)*log(2*l)+2*log(p1(x[1],l,m))+sum(log(p1(x[2:(N-1)],l,m)))+sum((S[1:N]-1)*log(l*p0m(trm.br[1:N],l,m))))  #I deleted the 2 in last log (20.2.2012)
 		if (survival==1){
 			res<- res + 2*log(1-m*p0m(x[1],l,m))}
 		}
 		res
    }
    if (lambda>0) {final<-dev(c(mu/lambda,lambda-mu)) } else {
    out <- subplex( c(0, 0.2),dev, hessian = TRUE)
    #res<-bd.groups.conf(out,dev)
    final<-out}
    final
    }
