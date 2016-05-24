manage <- function(init, epistep, vacgrid, costs, T=40, Tstop=T,
                   pinit = list(b=.1,k =.02,nu=.2,mu=.1),  ## initial parameter values
                   ## hyperparmeter values
                   hyper = list(bh=c(1,3), kh=c(1,3), nuh=c(1,1), muh=c(1,1)),
                   vac0 = list(frac=0, stop=0),  ## policy enacted before estimation, if any
                   MCvits=10, MCMCpits=1000, bkrate=1, vacsamps=100, ## MC sizes
                   start=8, ...)  ## time step to begin vaccination calculation
{
  ## get the initial last setting in epistep
  last <- formals(epistep)$last

  ## check for validity of costs argument
  if(is.null(costs) && !is.null(vacgrid))
    stop("costs cannot be NULL when vacgrid != NULL")

  ## check to ensure Tstop is <= than T
  if(Tstop > T) stop("must have final time horizon T >= Tstop")
  
  howmany <- 0
  doagain <- TRUE
  while(doagain && howmany<100)
    {      
      ## initialize the soln data frame
      soln <- init.soln(init, Tstop)
      
      ## initialize the prime data frame
      prime <- data.frame(matrix(0, ncol=2, nrow=Tstop))
      names(prime) <- c("S", "I")
      if(class(init) == "epiman") {
        ti <- nrow(init$prime)
        prime[1:ti,] <- init$prime
      } else { ti <- 1 }
      
      ## initialize the vaccination and culling policies
      VAC <- vac0$frac
      STOP <- vac0$stop
      allVfracs <- NULL  ## store the distn of vfracs given samples of parameter(s)
      allVstops <- NULL  ## similarly for vstops
      allVtimes <- NULL  ## and for vtimes
      allpolicies <- matrix(c(vac0$frac,vac0$stop), nrow=1)
      havesamps <- FALSE
      samp <- NULL
      ## starting values for the SIR parameters
      b <- pinit$b; k <- pinit$k; nu <- pinit$nu; mu <- pinit$mu;
      
      doagain <- FALSE
      alive <- 1
      ## step through time
      for( i in (ti+1):Tstop ){
        
        ## if you get to the start time, and the epidemic hasn't grown, start over
        if(i==start && soln$I[i-1]<=soln$I[1]) 
          {
            doagain <- TRUE
            break;
          }
        
        ## check if there is a need to vaccinate
        if(! (is.null(vacgrid) || soln$S[i-1] == 0 || soln$I[i-1] == 0) &&
           i>=start && havesamps)
          {
            ## call newvacpolicy so that vcsamps samples are taken
            ## conditional on a thinned version of samp
            VACs <- vaccinate(soln$S[i-1], soln$I[i-1], samp, vacsamps,
                               costs, MCvits, T-i+1, vacgrid)
            
            ## make the current vaccination policy equal to the mean
            mostbestpol <- which(VACs$bestC==max(VACs$bestC),arr.ind=TRUE)[1,]
            VAC <- vacgrid$fracs[mostbestpol[1]]
            STOP<- vacgrid$stops[mostbestpol[2]]

            if(soln$S[i-1]>STOP && VAC>0)
              {
                ## concatenate the policy for this time step to the history
                allVfracs <- rbind(allVfracs, VACs$Vfractions)
                allVstops <- rbind(allVstops, VACs$Vstoptimes)
                allVtimes <- c(allVtimes,i)
              }
            else
              {
                VAC<-0;STOP<-0;
              }
          }
        ## run epidemic one step forward
        allpolicies <- rbind(allpolicies, c(VAC,STOP))
        
        ## simulate one step forward in the epidemic, with epistep, with
        ## vaccinations and cullings coming first in epimanage
        out <- epimanage(soln=soln, epistep=epistep, i=i, VAC=VAC, STOP=STOP,
                         last=last, ...)
        last <- out$last
        out <- out$out
        
        ## update (prime) totals
        prime[i,] <- epi.update.prime(soln[i-1,], out)
        
        ## update (soln) totals
        soln[i,] <- epi.update.soln(soln[i-1,], out, costs)

        ## keep track of last time that the epidemic is "live"
        if((soln$S[i]>0 && soln$I[i]>0)) alive <- i
        
        ## do T MH/Gibbs steps to sample from the join distrib of b, k, mu and nu
        if(is.null(vacgrid) && i != Tstop) next;
        if(soln$S[i]>0 && soln$I[i]>0 || (i == Tstop)) {
          samp <- mcmc.bknumu(MCMCpits, bkrate, b, k, nu, mu, soln$itilde[2:alive],
                              soln$rtilde[2:alive], soln$dtilde[2:alive], prime$S[2:alive],
                              prime$I[2:alive], hyper)
          
          ## collect means    
          b <- mean(samp$b); k <- mean(samp$k); nu <- mean(samp$nu);
          mu <- mean(samp$mu)
          havesamps <- TRUE
        }
      }
    }	

  ## make allpolicies into a data frame
  allpolicies <- data.frame(allpolicies)
  names(allpolicies) <- c("frac", "stop")

  ## start to construct outputs
  r <- list(soln=soln, samp=samp, prime=prime, hyper=hyper)

  ## some outputs meaninglist when vacgrid = NULL
  if(is.null(vacgrid)) r$vachist <- r$vactimes <- r$pols <- NULL
  else {
    r$vachist <- list(fracs=allVfracs, stops=allVstops)
    r$vactimes <- allVtimes
    r$pols <- allpolicies
  }

  ## save call and return
  r$call <- match.call()
  class(r) <- "epiman"
  return(r)
}


## init.soln:
##
## function to initialize the solution -- the data frame
## that contains all of the information about the epidemic
## -- starting with a particular number of succeptables and
## infecteds, for a total time horizon t
init.soln <- function(init, T)
  {
    soln <- data.frame(matrix(0, ncol=11, nrow=T))
    names(soln) <- c("TIME", "S", "I", "R", "D", "itilde", "rtilde", "dtilde",
                   "V", "QC","C")

    if(class(init) == "epiman") {
      if(T < nrow(init$soln))
        stop("time horizon (T) shorter than rows in epiman object")
      soln[1:nrow(init$soln),] <- init$soln
    } else {
      soln$TIME[1] <- 1       # column 1: time
      soln$S[1] <- init$S     # column 2: susceptibles
      soln$I[1] <- init$I     # column 3: infecteds
    }

    return(soln)
  }


## epi.update.prime:
##
## the real number of succeptiables and infecteds that are
## around at the beginning of the current time step, where
## soln contains the state at the previous timestep and
## out contains information about the number just culled
## and just vaccinated
epi.update.prime <- function(soln, out)
  {
    prime <- data.frame(S=soln$S - out$justvacc,
                        I=soln$I - out$justcull)
    return(prime)
  }


## epi.update.soln:
##
## advance time by one unit in the soln data frame based
## on the previous time's state information (contained in s)
## and information (in out) about the number just culled,
## vaccinated, died, recovered, newly infected, etc
epi.update.soln <- function(s, out, costs)
  {
    if(!is.null(costs)) {
      cost.update <- s$C
      cost.update <- cost.update + costs$infect*(s$I-out$justdied)
      cost.update <- cost.update + costs$vac*out$justvacc
      cost.update <- cost.update + costs$death*out$justdied
    } else cost.update <- NA
    
    soln <-
      data.frame(
                 TIME = s$TIME+1,
                 S = s$S - out$justvacc - out$newi,
                 I = s$I - out$justcull - out$justrec - out$justdied + out$newi,
                 R = s$R + out$justrec, 
                 D = s$D + out$justdied,
                 itilde = out$newi,
                 rtilde = out$justrec,
                 dtilde = out$justdied,
                 V = s$V + out$justvacc,  
                 QC = s$QC + out$justcull,
                 C = cost.update)
    
    return(soln)
  }


## vaccinate:
##
## function to call newvacpolicy on a thinned version of the
## samples contained in samp.  The thinning is determined by
## the desired number of samples from the vaccination policy
## (vacsamps).  nI and nS are the number of infecteds and
## succeptibles in the current timestep.  Vold is the last
## used vaccination policy.  the rest of the arguments are
## identical to the manage function.
vaccinate <- function(nS,nI, samp, vacsamps, costs,
                      MCvits, timehorizon, vacgrid)
  {
    
    ## total number of active individuals
    N <- nI + nS

    ## indices of the thinned sample
    ll <- floor(seq(1,nrow(samp), length.out=vacsamps))

    ## the thinned sample
    use.samp <- samp[ll,]
    Vfractions <- rep(0,length(ll))
    Vstoptimes <- rep(0,length(ll))
    Vbestcosts <- matrix(0,length(vacgrid$fracs),length(vacgrid$stops))
    Vcosts <- array(0,c(length(ll),length(vacgrid$fracs),length(vacgrid$stops)))
    ## get the vaccination policy for each thinned sample,
    ## starting the next call where we left off on the
    ## last call
    for(j in 1:length(ll)){

      ## call the function which gets the sample
      VP <- VarStopTimePolicy(nS,nI,timehorizon,use.samp$b[j],
                              use.samp$k[j],use.samp$nu[j],use.samp$mu[j],
                              costs$vac, costs$death, costs$infect,
                              MCvits, vacgrid$fracs, vacgrid$stops,
                              midepidemic=TRUE,starttime=0)
      bestpol <- which(VP==min(VP),arr.ind=TRUE)[1,]
      Vfractions[j] <-vacgrid$fracs[bestpol[1]]
      Vstoptimes[j] <-vacgrid$stops[bestpol[2]]
      Vbestcosts[bestpol[1],bestpol[2]] <- Vbestcosts[bestpol[1],bestpol[2]]+1
      Vcosts[j,,] <- VP
      if(max(Vbestcosts)>vacsamps/2)
        {
          break
        }
    }

    return(list(Vfractions=Vfractions,Vstoptimes=Vstoptimes,C=Vcosts,bestC=Vbestcosts))
  }


## epistep:
##
## simulate the epidemic one time step foreward based on
## the previous time step and the parameters k, b, nu and mu
epistep <- function(SIR, last=NULL, true=list(b=0.00218, k=10, nu=0.4, mu=0))
  {
    b <- true$b; k <- true$k; nu <- true$nu; mu <- true$mu
    
    ## compute the number of removed
    rem <- rbinom(1, SIR$I, 1-exp(-(nu + mu)))
    
    ## compute the number of recovered
    rec <- rbinom(1, rem, nu/(nu + mu))
    
    ## compute the binomial probability of a new infection
    ## can use soln[i-1,3] or soln[i-1,3]-justcull
    p <- 1-(k/(k+b*(SIR$I)))^k   
    
    ## sample from the binomial distribution
    infect <- rbinom(1, SIR$S, p)
    
    ## compute the number which just died
    dead <- rem - rec
    
    return(list(rem=rem, rec=rec, infect=infect, dead=dead))
  }


## epimanage:
##
## simulate the epidemic one time step forward.  After managing
## via vaccinations and culling, then apply epistep on
## whatever is left.  The elipses argument (...) should contain
## whatever epistep needs to do its business, like b,k,nu, and mu
## (true) for example

epimanage <- function(soln, epistep=epistep, i, VAC=0.1, STOP=1, last=NULL, ...)
{ 
    ## compute the vacc/cull policy
    if(soln$S[i-1]>STOP) {
      vac <- ceiling(VAC*soln$S[i-1])
    } else vac <- 0
    
    cull<-0
    
    ## use epistep to simulate one time step forward based on the
    ## previous time step and the parameters b, k, nu, and mu
    SIR <- list(S=soln$S[i-1]-vac, I=soln$I[i-1]-cull, R=NA)
    just <- epistep(SIR, last=last, ...)

    ## assemble the outputs
    out <- matrix(0, ncol=5, nrow=1)
    out[1,1] <- vac
    out[1,2] <- cull
    out[1,3] <- just$rec
    out[1,4] <- just$infect
    out[1,5] <- just$dead
    out <- data.frame(out)
    names(out) <- c("justvacc", "justcull", "justrec", "newi", "justdied")

    return(list(out=out, last=just))
}
