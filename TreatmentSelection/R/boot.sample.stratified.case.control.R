boot.sample.stratified.case.control <-
function( event, trt,rho){

   n = length(trt)
   N.t0.d1 <- round(rho[1]*rho[3]) # how many cases with trt=0 in cohort
   N.t1.d1 <- round(rho[1]*rho[5]) # how many cases with trt=1 in cohort
   N.t1    <- round(rho[1]*(rho[4]+rho[5])) # how many with trt=1 in cohort
   N       <- rho[1]   
   f <- rho[6:7]

   N.t0 <- N - N.t1 #N is cohort sample size, so this is number of trt=0 in the cohort
   k.0 <- round(sum(event == 0 & trt == 0)/sum(event == 1 & trt == 0)) # number of controls per case in untreated arm
   k.1 <- round(sum(event == 0 & trt == 1)/sum(event == 1 & trt == 1)) #number of controls per case in treated arm
  
 #sample conditional on event and trt strata

   Nstar.trt1 <- rbinom(1, size = N, prob = N.t1/N) #sample number of treatments
   Nstar.trt0 <- N - Nstar.trt1 #number of controls for this boot

   Nstar.trt1.event1 <- rbinom(1, size = Nstar.trt1, prob = N.t1.d1/N.t1)  # sample number of controls in this treatment arm
   Nstar.trt0.event1 <- rbinom(1, size = Nstar.trt0, prob = N.t0.d1/N.t0)  # same 

   n.t0.d1 <- f[1]*Nstar.trt0.event1 # rewriting to homogenize notation
   n.t1.d1 <- f[2]*Nstar.trt1.event1

   n.t0.d0 <- k.0*n.t0.d1 #number of controls to sample in each trt arm
   n.t1.d0 <- k.1*n.t1.d1
    

    Resample <- TRUE
    tryN = 0
    while(Resample){
    
       myall   <- 1:n 
   boot.t0.d0 <- sample( myall[trt==0 & event==0], size = n.t0.d0,  replace = TRUE)
   boot.t1.d0 <- sample( myall[trt==1 & event==0], size = n.t1.d0,  replace = TRUE)
   boot.t0.d1 <- sample( myall[trt==0 & event==1], size = n.t0.d1,  replace = TRUE)
   boot.t1.d1 <- sample( myall[trt==1 & event==1], size = n.t1.d1,  replace = TRUE)
   ind <- c(boot.t0.d0, boot.t1.d0, boot.t0.d1, boot.t1.d1)
    
    Resample <- any(table(event[ind], trt[ind]) == 0) | length(unique(trt[ind]))==1 | length(unique(event[ind]))==1
    
    if(Resample){ warning("Bootstrap sample failed to sample individuals in each trt/event strata. Had to resample") ; tryN = tryN +1}
    if(tryN >10) stop("Had to resample too many bootstrap replicates due to too few observations in trt/event strata")
    }


   

   rho.b = c( N, (1 - Nstar.trt0.event1/Nstar.trt0)*Nstar.trt0/N,
                 Nstar.trt0.event1/N, 
                 (1 - Nstar.trt1.event1/Nstar.trt1)*Nstar.trt1/N,
                 Nstar.trt1.event1/N, f)
   #c((Nstar.trt0.event1)/(Nstar.trt0), (Nstar.trt1.event1)/(Nstar.trt1), Nstar.trt0.event1, Nstar.trt1.event1, Nstar.trt1, N, k)
   
   return(c(rho.b, ind))

}
