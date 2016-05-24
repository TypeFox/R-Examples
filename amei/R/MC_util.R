SimulateEpidemic <- function(init,T,params,vacc,vaccstop,costs,starttime, ...)
  {
    if(!is.function(params)) {
      soln = .C("RSimulateEpidemic",
        S = as.integer(rep(init$S,T)),
        I = as.integer(rep(init$I,T)),
        R = as.integer(rep(init$R,T)),
        D = as.integer(rep(init$D,T)),
        V = as.integer(rep(0,T)),
        C = as.double(rep(0,T)),
        T = as.integer(T),
        b = as.double(params$b),
        k = as.double(params$k),
        nu = as.double(params$nu),
        mu = as.double(params$mu),
        vacc = as.double(vacc),
        vaccstop = as.integer(vaccstop),
        cvacc = as.double(costs$vac),
        cdeath = as.double(costs$death),
        cinfected = as.double(costs$infect),
        starttime = as.integer(starttime),
        PACKAGE = "amei")
    } else {
      
      ## params is an epistep
      epistep <- params
      
      ## initialize the soln data frame
      soln <- init.soln(init, T)
      
      ## initialize the prime data frame
      prime <- data.frame(matrix(0, ncol=2, nrow=T))
      names(prime) <- c("S", "I")
      
      ## get the initial last setting in epistep
      last <- formals(epistep)$last
      
      for( i in 2:T ){
        
        ## deal with starttime
        if(i <= starttime) { VAC <- 0; STOP <- 0 }
        else { VAC <- vacc; STOP <- vaccstop }
        
        ## treat params like an epistep function
        out <- epimanage(soln=soln, epistep=epistep, i=i, VAC=VAC, STOP=STOP,
                         last=last, ...)
        last <- out$last
        out <- out$out
        
        ## update (prime) totals
        prime[i,] <- epi.update.prime(soln[i-1,], out)
        
        ## update (soln) totals
        soln[i,] <- epi.update.soln(soln[i-1,], out, costs)
      }

      
    }

    list(S = soln$S, I=soln$I, R=soln$R, D=soln$D, V=soln$V, C=soln$C)
  }

## ok, if midepidemic is true, then we're just going to accept all epidemics
## if midepidemic is false, then we're going to assume we're at the beginning 
## of an epidemic, and that starttime is something >0.  in this case, we'll 
## throw away any epidemics for which the number of infecteds hasn't grown 
## by startime. 
VarStopTimePolicy <- function(S0,I0,T,b,k,nu,mu,cvacc,cdeath,cinfected,
                             MCvits,Vprobs,Vstops,midepidemic,starttime)
  {
    cout = .C("RVarStopTimePolicy",
      S0 = as.integer(S0),
      I0 = as.integer(I0),
      T = as.integer(T),
      b = as.double(b),
      k = as.double(k),
      nu = as.double(nu),
      mu = as.double(mu),
      cvacc = as.double(cvacc),
      cdeath = as.double(cdeath),
      cinfected = as.double(cinfected),
      MCvits = as.integer(MCvits),
      Vprobs = as.double(Vprobs),
      nVprobs = as.integer(length(Vprobs)),
      Vstops = as.integer(Vstops),
      nVstops = as.integer(length(Vstops)),
      EC = double(length(Vprobs)*length(Vstops)),
      midepidemic = as.integer(midepidemic),
      starttime = as.integer(starttime),
      PACKAGE = "amei")
    C <- matrix(cout$EC,length(Vprobs),length(Vstops),byrow=TRUE)
    return(C)
  }


SimulateManagementQuantiles <- function(epistep,Time,init, pinit, hyper, vac0,
                                        costs, start, MCvits, MCMCpits, bkrate, vacsamps,
                                        vacgrid, nreps,lowerq, upperq, verb=FALSE, ...)
  {
    Sall <- matrix(0,nrow=Time,ncol=nreps)
    Iall <- matrix(0,nrow=Time,ncol=nreps)
    Rall <- matrix(0,nrow=Time,ncol=nreps)
    Dall <- matrix(0,nrow=Time,ncol=nreps)
    Vall <- matrix(0,nrow=Time,ncol=nreps)
    Call <- matrix(0,Time,nreps)
    PoliciesAll <- array(0,c(Time,2,nreps))
    
    for(n in 1:nreps)
      {
        if(verb) cat("*** Simulating epidemic",n,"***\n")
        foo <- manage(epistep=epistep,pinit=pinit,T=Time,Tstop=Time,init=init,
                      hyper=hyper, vac0, costs=costs, MCMCpits=MCMCpits, bkrate=bkrate,
                      vacsamps=vacsamps, vacgrid=vacgrid, start=start, ...)
        Sall[,n] <- foo$soln$S
        Iall[,n] <- foo$soln$I
        Rall[,n] <- foo$soln$R
        Dall[,n] <- foo$soln$D
        Vall[,n] <- foo$soln$V
        Call[,n] <- foo$soln$C
        PoliciesAll[,,n] <- as.matrix(foo$pols)
      }
    
    SQ1 <- apply(Sall,1,quantile,prob=lowerq)
    Smean <- apply(Sall,1,mean)
    Smed <- apply(Sall,1,median)
    SQ3 <- apply(Sall,1,quantile,prob=upperq)
    
    IQ1 <- apply(Iall,1,quantile,prob=lowerq)
    Imean <- apply(Iall,1,mean)
    Imed <- apply(Iall,1,median)
    IQ3 <- apply(Iall,1,quantile,prob=upperq)
    
    RQ1 <- apply(Rall,1,quantile,prob=lowerq)
    Rmean <- apply(Rall,1,mean)
    Rmed <- apply(Rall,1,median)
    RQ3 <- apply(Rall,1,quantile,prob=upperq)
    
    DQ1 <- apply(Dall,1,quantile,prob=lowerq)
    Dmean <- apply(Dall,1,mean)
    Dmed <- apply(Dall,1,median)
    DQ3 <- apply(Dall,1,quantile,prob=upperq)
    
    VQ1 <- apply(Vall,1,quantile,prob=lowerq)
    Vmean <- apply(Vall,1,mean)
    Vmed <- apply(Vall,1,median)
    VQ3 <- apply(Vall,1,quantile,prob=upperq)
    
    CQ1 <- apply(Call,1,quantile,prob=lowerq)
    Cmean <- apply(Call,1,mean)
    Cmed <- apply(Call,1,median)
    CQ3 <- apply(Call,1,quantile,prob=upperq)
    
    PolQ1 <- apply(PoliciesAll,c(1,2),quantile,prob=lowerq)
    Polmean <- apply(PoliciesAll,c(1,2),mean)
    Polmed <- apply(PoliciesAll,c(1,2),median)
    PolQ3 <- apply(PoliciesAll,c(1,2),quantile,prob=upperq)
    
    list(Q1 = data.frame(S=SQ1,I=IQ1,R=RQ1,D=DQ1,V=VQ1,C=CQ1,frac=PolQ1[,1],stop=PolQ1[,2]),
         Mean = data.frame(S=Smean,I=Imean,R=Rmean,D=Dmean,V=Vmean,C=Cmean,frac=Polmean[,1],stop=Polmean[,2]),
         Median = data.frame(S=Smed,I=Imed,R=Rmed,D=Dmed,V=Vmed,C=Cmed,frac=Polmed[,1],stop=Polmed[,2]),
         Q3 = data.frame(S=SQ3,I=IQ3,R=RQ3,D=DQ3,V=VQ3,C=CQ3,frac=PolQ3[,1],stop=PolQ3[,2]))       
  }

SimulateEpidemicQuantiles <- function(init,T,params,vacc,vaccstop,costs,
  nreps,lowerq,upperq,midepidemic,starttime)
  {
    Sall <- matrix(0,nrow=T,ncol=nreps)
    Iall <- matrix(0,nrow=T,ncol=nreps)
    Rall <- matrix(0,nrow=T,ncol=nreps)
    Dall <- matrix(0,nrow=T,ncol=nreps)
    Vall <- matrix(0,nrow=T,ncol=nreps)
    Call <- matrix(0,T,nreps)
    for(n in 1:nreps)
      {
        isvalid <- TRUE;isvalidcount<-0
        while(isvalid)
          {
            tmpsim <- SimulateEpidemic(init,T,params,vacc,vaccstop,costs,starttime)
            if(!midepidemic)
              {
                isvalidcount <- isvalidcount+1
                if(tmpsim$I[starttime-1]>init$I)
                  isvalid <- FALSE
                if(isvalidcount==100)
                  {
                    cat("Warning: <1% chance of an epidemic\n")
                    isvalid<-FALSE
                  }
              }
                
          }
        Sall[,n] <- tmpsim$S
        Iall[,n] <- tmpsim$I
        Rall[,n] <- tmpsim$R
        Dall[,n] <- tmpsim$D
        Vall[,n] <- tmpsim$V
        Call[,n] <- tmpsim$C
      }

    SQ1 <- apply(Sall,1,quantile,prob=lowerq)
    Smean <- apply(Sall,1,mean)
    Smed <- apply(Sall,1,median)
    SQ3 <- apply(Sall,1,quantile,prob=upperq)

    IQ1 <- apply(Iall,1,quantile,prob=lowerq)
    Imean <- apply(Iall,1,mean)
    Imed <- apply(Iall,1,median)
    IQ3 <- apply(Iall,1,quantile,prob=upperq)

    RQ1 <- apply(Rall,1,quantile,prob=lowerq)
    Rmean <- apply(Rall,1,mean)
    Rmed <- apply(Rall,1,median)
    RQ3 <- apply(Rall,1,quantile,prob=upperq)

    DQ1 <- apply(Dall,1,quantile,prob=lowerq)
    Dmean <- apply(Dall,1,mean)
    Dmed <- apply(Dall,1,median)
    DQ3 <- apply(Dall,1,quantile,prob=upperq)
    
    VQ1 <- apply(Vall,1,quantile,prob=lowerq)
    Vmean <- apply(Vall,1,mean)
    Vmed <- apply(Vall,1,median)
    VQ3 <- apply(Vall,1,quantile,prob=upperq)

    CQ1 <- apply(Call,1,quantile,prob=lowerq)
    Cmean <- apply(Call,1,mean)
    Cmed <- apply(Call,1,median)
    CQ3 <- apply(Call,1,quantile,prob=upperq)

    
    list(Q1 = data.frame(S=SQ1,I=IQ1,R=RQ1,D=DQ1,V=VQ1,C=CQ1),
         Mean = data.frame(S=Smean,I=Imean,R=Rmean,D=Dmean,V=Vmean,C=Cmean),
         Median = data.frame(S=Smed,I=Imed,R=Rmed,D=Dmed,V=Vmed,C=Cmed),
         Q3 = data.frame(S=SQ3,I=IQ3,R=RQ3,D=DQ3,V=VQ3,C=CQ3))       
  }
