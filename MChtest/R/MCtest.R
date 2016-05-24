"MCtest" <-
function(x,statistic,resample,bound=MCbound.precalc1,extreme="geq",seed=1234325){
    set.seed(seed)
    if (bound$type=="fixed"){ 
        out<-MCtest.fixed(x,statistic,resample,Nmax=bound$parms[1],extreme=extreme,conf.level=bound$conf.level,seed=seed) 
    }
    else {        
        out<-MCtest.default(x,statistic,resample,bound=bound,extreme=extreme) 
    }
    out
}

