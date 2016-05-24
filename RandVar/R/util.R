.getImageDistr <- function(f, distr){ 
    if (is(distr, "DiscreteDistribution"))
        return(DiscreteDistribution(prob=d(distr)(support(distr)), supp=f(support(distr))))
    if(is(try(return(f(distr)), silent = TRUE), "try-error")){
        rl <- function(n){ 
            xr <- r(distr)(n) 
            f(xr) 
        }

        n <- 10^getdistrOption("RtoDPQ.e")+1
        u <- seq(0,1,length=n+1); u <- (u[1:n]+u[2:(n+1)])/2
        qd <- q(distr)
        y <- f(qd(u))
    
        wmdn <- getdistrOption("warn.makeDNew")
        on.exit(distroptions(warn.makeDNew=wmdn))
        distroptions(warn.makeDNew=FALSE)
        
        if(length(unique(c(rl(10000),y)))==10000+length(y)){
           DPQnew <- RtoDPQ(r=rl, y=y)
           return(AbscontDistribution(r = rl, d = DPQnew$d, p = DPQnew$p, 
                                      q = DPQnew$q, .withArith = TRUE, 
                                      .withSim = TRUE, withgaps = FALSE))
        
        }else
           return(UnivarLebDecDistribution(r = rl, y = y))
    }
}
