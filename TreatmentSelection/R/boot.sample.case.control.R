boot.sample.case.control <-
function( event, trt, rho ){
    # f is the fraction of cases we sampled 
    #sample conditional on event
    f <- rho[4]
    N <- rho[1] #rho = c(Pr(event=1), cohort sample size, Pr(event==1|trt=0), Pr(event=1|trt=1), nmatch, No of trt==1)
    k <- round(sum(event==0)/sum(event==1)) #number of controls per case
   # N.t <- round(rho[2]*rho[1])

    N.d1.star <- rbinom(1, size = N, prob = (rho[3]))
    
    n.d1 <- f*N.d1.star #how many cases to sample

    n.d0 <- k*n.d1 #how many controls to sample. k is number of controls/cases we want

    
    myall   <- 1:length(trt) 
    
    Resample <- TRUE
    tryN = 0
    while(Resample){
    ind.d0 <- sample(myall[event==0], size = n.d0, replace = TRUE)
    ind.d1 <- sample(myall[event==1], size = n.d1, replace = TRUE)

    ind <- c(ind.d0, ind.d1)
        Resample <- any(table(event[ind], trt[ind]) == 0) | length(unique(trt[ind]))==1 | length(unique(event[ind]))==1
        if(Resample){ warning("Bootstrap sample failed to sample individuals in each trt/event strata. Had to resample") ; tryN = tryN +1}
    if(tryN >10) stop("Had to resample too many bootstrap replicates due to too few observations in trt/event strata")
    }
    
   # event.b <- event[ind]
   # trt.b <- trt[ind]

   # Pr.D1.cond.trt1.b <- mean(event.b[trt.b==1])*sum(trt.b==1)/N.t
   # Pr.D1.cond.trt0.b <- mean(event.b[trt.b==0])*sum(trt.b==0)/(N-N.t)
    

    rho.b <- c( N, rho[2], (N.d1.star)/N, f, -9999, -9999, -9999)  
    
    return(c(rho.b, ind))
}
