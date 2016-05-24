inclusionprobastrata<-function (strata, nh) 
{   N = length(strata)
    EPS = 1e-6
    if (min(unique(strata)) < 1) 
        stop("the stratification variable has incorect values (less than 1)\n")
    Nh=as.vector(table(strata))
    if(any(nh/Nh>1+EPS)) warning("in a stratum the sample size is larger than the population size\n")
    pik=nh[strata]/Nh[strata]
    pik
}
