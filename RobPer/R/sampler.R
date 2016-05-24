sampler <- function(ttype, npoints, ncycles, ps=1) {
	
    if(ttype=="equi") {
        tt <- seq(from=0, to =ps*ncycles, length.out=npoints+1)[-1]
    } else {
        unif.phase <- runif(npoints)
    }
    # Phase in [0,1]
	
    if(ttype=="sine") {
        cdf <- function(x) return(x-(cos(2*pi*x)-1)/(2*pi))
        cdf.norm <- function(x) {
            f <- cdf(x)-unif.phase
            f
        }
        phase <- BBsolve(par=rep(ps/2, npoints), fn=cdf.norm)$par
    }
    if(ttype=="trian") {
        phase <- sqrt(2/3*unif.phase)*(unif.phase<=2/3)-(sqrt((1-unif.phase)/3)-1)*(unif.phase>2/3)
    }
	
    # Cycle in {1,...,ncycles}
    if(ttype=="unif") {
        tt <- unif.phase*ncycles*ps
    }
    if(ttype%in%c("sine", "trian")) {
        random.cycles <- sample(ncycles, npoints, T)
        tt <- ps*(random.cycles-1+phase)
    }
    return(sort(tt))
}
