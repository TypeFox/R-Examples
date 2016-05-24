library(lgcp)

n <- 100
accs <- rep(0.5,n)

mcmcloop <- mcmcLoop(N=n,burnin=floor(n/10),thin=1,progressor=mcmcProgressTextBar) # it is important that this object is named "mcmcloop"

adsch <- andrieuthomsh(inith=1,alpha=1/2,C=1,targetacceptance=0.574)

h <- initialiseAMCMC(adsch)
hrec <- h
while(nextStep(mcmcloop)){
    ac <- accs[iteration(mcmcloop)] # it is important that this objects is named "ac"
    h <- updateAMCMC(adsch)
    hrec <- c(hrec,h)
}	

halt <- 1
haltrec <- halt
itno <- 1
while(itno<=100){
    itr <- itno
    if(itno>floor(n/10)){
        itr <- itno - floor(n/10)
    }
    halt <- exp(log(halt) + (1/(itr^(1/2)))*(accs[itr]-0.574))
    haltrec <- c(haltrec,halt)
    itno <- itno + 1
}

print(hrec)
print("")
print(haltrec)
print(all(hrec==haltrec))