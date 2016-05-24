##"_helpers" don't get exported in package namespace

                                        #assume gamma prior on lambda, mu, with parameters (alpha.L, beta.L) and
                                        # (alpha.M, beta.M) respectively.  the betas are "rates" (not "means").
                                        #Right now, nu = 0.  Would have to play with prior to get something else to work
                                        #data should be a partially observed CTMC (observed at discrete time points)
                                        # with times, states, and T as objects.  (times[1]=0).
                                        #N is number of full iterations to run the MCMC for.
                                        #burnIn is number to be discarded at beginning.
BD.MCMC.SC.CTMC_PO_1 <- function(Lguess, Mguess,
                                 beta.immig, #not a guess, of course
                                 alpha.L, beta.L, alpha.M, beta.M,
                                 data,
                                 burnIn=100, N=1000, n.fft=1024,
                                 verbose=1, verbFile=NULL,simMethod=-1,...){
  if (!is.null(verbFile)) sink(verbFile)
  if ( inherits(data, "CTMC_PO_1")){
    T <- getT(data);
  }
  else
    T <- data$T;##just so that can have identical function for ctmcpo1 and list
  param.history <- matrix(data=NA, nrow=N+1, ncol=2);
  ##fullBD.history <- vector(mode="list") ;#right now am not keeping the augmented states around
  param.history[1,] <- c(Lguess, Mguess);
  L <- Lguess;
  M <- Mguess;
  if (simMethod==0){
    simfn <- function(la,mu){bdARsimCondEnd(Naccepted=1, Ntotal=NULL, bd.PO=data,L=la, m=mu, nu=beta.immig*la)[[1]];}
  }
  else if (simMethod==1){
    simfn <- function(la,mu){sim.condBD(1,bd.PO=data, L=la, m=mu, nu=beta.immig*la, n.fft=n.fft)[[1]];}
  }
  else ##preferred
    simfn <- function(la,mu){chooseSim.condBD.1.CTMC_PO_1(bd.PO=data, L=la, m=mu, nu=beta.immig*la, n.fft=n.fft);}
  for (i in 1:N){
    BD.full <- simfn(L,M);
    ##BD.full <-  sim.condBD(1,bd.PO=data, L=L, m=M, nu=beta.immig*L, n.fft=n.fft)[[1]];
    ##BD.full <- bdARsimCondEnd(Naccepted=1, Ntotal=NULL, bd.PO=data,
    ##                       L=L, m=M, nu=beta.immig*L)[[1]];
    fullSummStats <- BDsummaryStats(BD.full);
    L <- rgamma(n=1, shape=alpha.L + fullSummStats["Nplus"],
                rate= beta.L + fullSummStats["Holdtime"] + beta.immig*T);
    M <- rgamma(n=1, shape=alpha.M + fullSummStats["Nminus"],
                rate = beta.M + fullSummStats["Holdtime"]);
    param.history[i+1,] <- c(L, M);
    if (verbose>0 && (i %% 30==0)){
      print(paste("BD.MCMC.SC: On the ",i,"th iteration params are ",
                  L,M,sep=" "));
    }
  }
  if (!is.null(verbFile)) sink(NULL);
  param.history[(burnIn+1):(N+1),];
}

BD.MCMC.SC.list <- BD.MCMC.SC.CTMC_PO_1;

BD.MCMC.SC.CTMC_PO_many <- function(Lguess, Mguess,
                                    beta.immig, #not a guess, of course
                                    alpha.L, beta.L, alpha.M, beta.M,
                                    data,
                                    burnIn=100, N=1000, n.fft=1024,
                                    verbose=1, verbFile=NULL,
                                    simMethod=-1,...){
  if (!is.null(verbFile)) sink(verbFile)
  T <- getT(data);
  param.history <- matrix(data=NA, nrow=N+1, ncol=2);
  ##fullBD.history <- vector(mode="list") ;#right now am not keeping the augmented states around
  param.history[1,] <- c(Lguess, Mguess);
  L <- Lguess;
  M <- Mguess;
  BD.fulls <- vector("list", length=length(data@BDMCsPO));
  if (simMethod==0){
    simfn <- function(ctmc1, la,mu){bdARsimCondEnd(Naccepted=1, Ntotal=NULL, bd.PO=ctmc1,L=la, m=mu, nu=beta.immig*la)[[1]];}
  }
  else if (simMethod==1){
    simfn <- function(ctmc1, la,mu){sim.condBD(1,bd.PO=ctmc1, L=la, m=mu, nu=beta.immig*la, n.fft=n.fft)[[1]];}
  }
  else ##PREFERRED:
    simfn <- function(ctmc1, la,mu){chooseSim.condBD.1.CTMC_PO_1(bd.PO=ctmc1, L=la, m=mu, nu=beta.immig*la, n.fft=n.fft);}
  for (i in 1:N){
    mySimmer <- function(ctmc1){
      ##sim.condBD(1,bd.PO=ctmc1, L=L, m=M, nu=beta.immig*L, n.fft=n.fft)[[1]]
      simfn(ctmc1=ctmc1,L,M)
    }
##########replace 'tryCatch' version with regular one 
    ##BD.fulls <- lapply(data@BDMCsPO, mySimmer);
        for (j in 1:length(data@BDMCsPO)){
          BD.fulls[[j]] <- tryCatch(mySimmer(data@BDMCsPO[[j]]),
                               error=function(e){
                                 print(paste("j is",j));
                                 print("The offending ctmc-po-1 entry of the data vector is")
                                 print(data@BDMCsPO[[j]]);
                                 e;
                               },
                               finally=NULL)
          ####BD.fulls[[j]] <- mySimmer(data@BDMCsPO[[j]]);
        }
    fullSummStatsMat <- sapply(BD.fulls, BDsummaryStats); ##nx3 matrix
    fullSummStats <- apply(fullSummStatsMat, 1, sum)
    L.old <- L; M.old <- M; ##entirely for debugging.
    L <- rgamma(n=1, shape=alpha.L + fullSummStats["Nplus"],
                rate= beta.L + fullSummStats["Holdtime"] + beta.immig*T);
    M <- rgamma(n=1, shape=alpha.M + fullSummStats["Nminus"],
                rate = beta.M + fullSummStats["Holdtime"]);

    ###DEBUG HERE
    if (L< 1e-10 || M< 1e-10){
      print("L or M is smaller than expected. here are params.")
      print(paste(alpha.L, beta.L,alpha.M,beta.M))
      print(paste(i,L,M,beta.immig,T))
      print(paste(L.old,M.old));
      print(fullSummStats)
      
    }


    ##     ######### HERE HERE DEBUG CODE
    ##     if (L < 1e-11 || M < 1e-11){
    ##       sink("BDMCMCSCCTMCPOmany",append=TRUE)
    ##       print("BD.MCMC.SC.CTMC_PO_many: L or M is confusingly small.")
    ##       print(paste("We are in iteration", i))
    ##       print(paste("L,M are", L,M));
    ##       print(paste("alpha.L,Nplus,betaL,Holdtime,beta,T,alpha.M,Nminus,betaM,Holdtime are",
    ##                   alpha.L, fullSummStats["Nplus"], beta.L, fullSummStats["Holdtime"],
    ##                   beta.immig,T,alpha.M,fullSummStats["Nminus"],
    ##                   beta.M,fullSummStats["Holdtime"]))
    ##       print(BD.fulls);
    ##       sink(NULL)
    ##       save(RNGvec, file="BDMCMCSCCTMCPOmany.RNGvec.rsav")
    ##       stop("quitting b/c got too small of param value");
    ##     }
    ##     ######## end HERE HERE DEBUG
    
    param.history[i+1,] <- c(L, M);
    if (verbose>0 && (i %% 30==0)){
      print(paste("BD.MCMC.SC: On the ",i,"th iteration params are ",
                  L,M,sep=" "));
    }
  }
  if (!is.null(verbFile)) sink(NULL)
  param.history[(burnIn+1):(N+1),];
}




####################
####################  Code for full BDI model MCMC
####################  Not exported because too slow:
####################   Never figured out how to quickly calculate coefficients.







