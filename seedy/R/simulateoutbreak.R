simulateoutbreak <-
function(init.sus, inf.rate, rem.rate, mut.rate, nmat=NULL, equi.pop=10000, shape=flat,
         init.inf=1, inoc.size=1, samples.per.time=1, samp.schedule="random", 
         samp.freq=500, full=FALSE, mincases=1, feedback=500, glen=100000, 
         ref.strain=NULL, ...) {
  
  # WARNINGS
  
  if (init.sus%%1!=0 || init.sus<1) {
    stop("Initial number of susceptibles must be a postive integer")
  }
  if (init.inf%%1!=0 || init.inf<1) {
    stop("Initial number of susceptibles must be a postive integer")
  }
  if (inoc.size%%1!=0 || inoc.size<1) {
    stop("Inoculum size must be a postive integer")
  }
  if (samp.freq%%1!=0 || samp.freq<1) {
    stop("samp.freq must be a postive integer")
  }
  if (feedback%%1!=0 || feedback<1) {
    stop("feedback must be a postive integer")
  }
  if (glen%%1!=0 || glen<1) {
    stop("Genome length must be a postive integer")
  }
  if (mincases%%1!=0 || mincases<1 || mincases>init.sus+init.inf) {
    stop("Minimum cases must be a postive integer not greater than init.sus+init.inf")
  }
  if (samples.per.time%%1!=0 || samples.per.time<1) {
    stop("samples.per.time must be a postive integer")
  }
  if (inf.rate<=0) {
    stop("Infection rate must be greater than zero")
  }
  if (rem.rate<=0) {
    stop("Removal rate must be greater than zero")
  }
  if (mut.rate<0 || mut.rate>=1) {
    stop("Mutation rate must be between 0 and 1")
  }
  if (!is.function(shape)) {
    stop("'shape' must be a function")
  } 
  if (!is.null(nmat)) {
    if (!is.matrix(nmat)) {
      stop("nmat must be a matrix")
    } else {
      if (nrow(nmat)!=init.sus+init.inf || ncol(nmat)!=init.sus+init.inf) {
        stop("nmat must have init.sus+init.inf rows and columns")
      } else if (sum(nmat<0)>0) {
        stop("All entries of nmat must be >= 0")
      }
    }
  }
  if (equi.pop%%1!=0 || equi.pop<=0) {
    stop("Equilibrium population size must be a postive integer")
  }
  if (!samp.schedule%in%c("random", "calendar", "individual")) {
    stop("samp.schedule must be 'random', 'calendar', or 'individual'")
  }
  
  #########################
  
  cat("\nSimulating outbreak:\n")
  cat("N=",equi.pop, ", b=", inoc.size, ", beta=", inf.rate, 
      ", gamma=", rem.rate, ", mu=", mut.rate, "\n\n", sep="")
  trigger <- FALSE
  
  at <- 0
  while (!trigger) { # Repeat if < mincases are infected
    newinfect <- 0
    at <- at+1
    cat("Attempt ", at, "\n", sep="")
    eff.cur.inf <- NULL
    cur.inf <- 1:init.inf # vector of infected person IDs
    cur.sus <- init.inf+(1:init.sus) # vector of susceptible person IDs
    time <- 0 # in bacterial generations
    inf.times <- rep(0,init.inf) # vector of infection times
    rec.times <- rgeom(init.inf,rem.rate)+2 # vector of removal times
    tot.inf <- init.inf
    inf.ID <- 1:init.inf
    if (samp.schedule=="random") {
      sample.times <- NULL
      for (i in 1:init.inf) {
        sample.times <- c(sample.times, sample(1:(rec.times[i]-1),1)) # vector of sampling times
      }
    } else {
      sample.times <- rep(samp.freq, init.inf)
    }
    inf.source <- rep(0,init.inf) # source of infection for each individual
    if (is.null(ref.strain)) {
      ref.strain <- sample(1:4, glen, replace=T) # reference strain
    } else {
      glen <- length(ref.strain)
    }
    totcurstrains <- 1 # current list of strains
    uniquestrains <- 1 # Number of unique strain types
    
    libr <- list() # list of mutation locations for each genotype
    mut.nuc <- list() # nucleotides at mutation locations
    freq.log <- list() # List of strain frequencies for each infective
    strain.log <- list() # Strain IDs for each within host population
    
    # Initialize logs
    for (i in 1:init.inf) {
      if (sample.times[i]>rec.times[i]) {
        sample.times[i] <- Inf
      }
      if (i == 1) {
        libr[[i]] <- NA
        mut.nuc[[i]] <- NA      
      } else {
        libr[[i]] <- sample(glen,1)
        mut.nuc[[i]] <- sample((1:4)[-ref.strain[libr[[i]]]], 1)
      }
      freq.log[[i]] <- 1
      strain.log[[i]] <- 1
    }
    
    for (i in (init.inf+1):(init.sus+init.inf)) {
      freq.log[[i]] <- 0
      strain.log[[i]] <- 0
    }
    
    current.infected <- init.inf
    types <- 1 # Cumulative number of strain types
    
    #Sample logs
    if (full) {
      obs.freq <- list()
      obs.strain <- list()
      pID <- NULL
    } else {
      sampleWGS <- NULL
      samplepick <- NULL
    }
    
    sampletimes <- NULL
    sampleID <- NULL
    
    while (length(cur.inf) > 0) { # Cycle through bacterial generations until epidemic ceases
      time <- time+1
      if (time%in%rec.times) { # recovery?
        recover <- inf.ID[which(rec.times==time)] # who has recovered?
        cur.inf <- cur.inf[-which(cur.inf%in%recover)] # remove infective(s)
        for (r in 1:length(recover)) {
          strain.log[[recover[r]]] <- 0
          freq.log[[recover[r]]] <- 0
        }
        if (length(cur.inf)==0) { # If no more infectives
          cat("t=", time, ", S=", length(cur.sus), ", I=", length(cur.inf), 
              ", total genotypes=0\n", sep="")
        }
      }
      if (time%%feedback==0 && length(cur.inf)>0) { # output current status every x generations
        cat("t=", time, ", S=", length(cur.sus), ", I=", length(cur.inf), 
            ", total genotypes=", length(unique(as.numeric(unlist(strain.log)))), ", next rec time=", 
            min(rec.times[which(rec.times>time)]), sep="")
        if (length(cur.sus)>0) {
          cat("\n")
        } else {
          cat(", final removal time=", max(rec.times), "\n", sep="")
        }
      }
      # calculate force of infection
      
      if (is.null(nmat)) {
        curinfrate <- inf.rate*length(cur.inf)/(init.inf+init.sus)
        pinf <- rbinom(length(cur.sus),1,curinfrate)
      } else {
        inf.mat <- nmat*inf.rate/(init.inf+init.sus)
        if (length(cur.inf)>1) {
        	curinfrate <- apply(inf.mat[cur.sus,cur.inf],1,sum)
        } else {
        	curinfrate <- inf.mat[cur.sus,]
        }
        pinf <- rbinom(length(cur.sus),1,curinfrate)
      }
      
      #if (runif(1,0,1) < curinfrate) { # infection?
      if (sum(pinf)>0) {
      	for (kp in 1:sum(pinf)) {
	        tot.inf <- tot.inf+1
	        # who got infected? and by whom?
	        newinfect <- cur.sus[which(pinf==1)[kp]]
	          
	        if (length(cur.inf)==1) {
	            inf.source <- c(inf.source, cur.inf)
	        } else if (is.null(nmat)) {
	            inf.source <- c(inf.source, sample(cur.inf,1)) # sample source at random
	        } else {
	          	inf.source <- c(inf.source, sample(cur.inf, 1, prob=inf.mat[newinfect,cur.inf]))
	        }
	        
	        inf.ID <- c(inf.ID, newinfect)
	        cur.inf <- c(cur.inf, newinfect) # add to current infectives
	        cur.sus <- cur.sus[-which(cur.sus==newinfect)] # Remove susceptible
	        inf.times <- c(inf.times, time) # new infection time
	        rec.times <- c(rec.times, time+max(1,rgeom(1,rem.rate))) # Don't recover today!
	        if (samp.schedule == "individual") {
	          sample.times <- c(sample.times, time+samp.freq)
	        } else if (samp.schedule == "calendar") {
	          sample.times <- c(sample.times, ceiling(time/samp.freq)*samp.freq)
	        } else if (samp.schedule == "random") {
	          if (rec.times[tot.inf]>time+1) {
	            sample.times <- c(sample.times, sample(time:(rec.times[tot.inf]-1),1))
	          } else {
	            sample.times <- c(sample.times, time)
	          }
	        }
	        if (sample.times[tot.inf]>=rec.times[tot.inf]) {
	          sample.times[tot.inf] <- Inf
	        }
	        # pass on strain
	        src <- inf.ID[which(inf.ID==inf.source[tot.inf])] # Source of infection
	        if (length(strain.log[[src]])==1) { # if source has clonal infection
	          inoc.samp <- rep(strain.log[[src]], inoc.size)
	          if (0%in%inoc.samp) {
	            stop("Zeroes in inoculum")
	          }
	        } else {
	          inoc.samp <- sample(strain.log[[src]], inoc.size, 
	                              prob=freq.log[[src]], replace=T) # take random sample
	          if (0%in%inoc.samp) {
	            stop("Zeroes in inoculum")
	          }
	        }
	        strain.log[[newinfect]] <- unique(inoc.samp) # distinct types in new infection
	        f <- numeric(length(unique(inoc.samp)))
	        k <- 1
	        for (i in unique(inoc.samp)) {
	          f[k] <- sum(inoc.samp==i)
	          k <- k+1
	        }
	        freq.log[[newinfect]] <- f # frequency of types
      	}
      }
      
      # mutate existing strains for each individual
      if (is.null(eff.cur.inf)) {
        cinf <- inf.ID[which(inf.ID%in%cur.inf)]
      } else {
        cinf <- eff.cur.inf
      }
      for (i in cinf) {
        pop.size <- sum(freq.log[[i]])
        death.prob <- min(0.5 + 0.5*(pop.size-shape(time,span=rec.times[i]-inf.times[i]+1,equi.pop,...))/shape(time,span=rec.times[i]-inf.times[i]+1,equi.pop,...),1)
        if (length(freq.log)==0 || sum(is.na(freq.log))>0 || death.prob>1 || death.prob<0) {
          cat("deathprob=", death.prob, "\npop.size=", pop.size, "\nequi.pop=", equi.pop, "\nFreq.log:\n")
          if (length(freq.log)>0) {
            for (k in 1:length(freq.log)) {
              cat(freq.log[[i]][k], "\n")
            }
          }
        }
        freq.log[[i]] <- 2*rbinom(length(freq.log[[i]]), freq.log[[i]], 1-death.prob)
        if (0 %in% freq.log[[i]]) {
          zeros <- which(freq.log[[i]]==0)
          if (length(zeros)==length(freq.log[[i]])) {
            freq.log[[i]] <- 1
            strain.log[[i]] <- strain.log[[i]][1]
          } else {
            freq.log[[i]] <- freq.log[[i]][-zeros]
            strain.log[[i]] <- strain.log[[i]][-zeros]
          }
        }
        if (length(strain.log[[i]])!=length(freq.log[[i]])) {
          stop("Error")
        }
        n.mutations <- rbinom(1, sum(freq.log[[i]]), mut.rate)
        if (n.mutations > 0) {
          for (mt in 1:n.mutations) {
            types <- types+1
            if (length(strain.log[[i]])==1) {
              mutate.grp <- strain.log[[i]]
            } else {
              mutate.grp <- sample(strain.log[[i]], 1, prob=freq.log[[i]])
            }
            if (mutate.grp %in% totcurstrains) {
              mut.loc <- sample(glen, 1)
              mut.nuc[[length(totcurstrains)+1]] <- 
                mut.nuc[[which(totcurstrains==mutate.grp)]]
              if (mut.loc %in% libr[[which(totcurstrains==mutate.grp)]]) { # if mutation at existing location
                kn <- which(libr[[which(totcurstrains==mutate.grp)]]==mut.loc)
                mut.nuc[[length(totcurstrains)+1]][kn] <- sample((1:4)[-mut.nuc[[which(totcurstrains==mutate.grp)]][kn]], 1)
                libr[[length(totcurstrains)+1]] <- libr[[which(totcurstrains==mutate.grp)]]
              } else {
                mut.nuc[[length(totcurstrains)+1]] <- 
                  c(mut.nuc[[length(totcurstrains)+1]], 
                    sample((1:4)[-ref.strain[mut.loc]], 1))
                libr[[length(totcurstrains)+1]] <- 
                  c(libr[[which(totcurstrains==mutate.grp)]], mut.loc)
              }
              if (sum(is.na(mut.nuc[[length(totcurstrains)+1]]))>0) {
                mut.nuc[[length(totcurstrains)+1]] <- mut.nuc[[length(totcurstrains)+1]][-is.na(mut.nuc[[length(totcurstrains)+1]])]
                libr[[length(totcurstrains)+1]] <- libr[[length(totcurstrains)+1]][-is.na(libr[[length(totcurstrains)+1]])]
              }
              strain.log[[i]] <- c(strain.log[[i]], types)
              freq.log[[i]] <- c(freq.log[[i]], 1)
              totcurstrains <- c(totcurstrains, types)
            }
          }
        }
      }
      # take samples, make observations
      if (time%in%sample.times) {
        smpat <- inf.ID[which(sample.times==time)]
        for (i in smpat) {
          if (full) {
            n <- length(obs.freq)+1
            obs.freq[[n]] <- freq.log[[i]]
            obs.strain[[n]] <- strain.log[[i]]
            sampleID <- c(sampleID, n)
            pID <- c(pID, i)
            sampletimes <- c(sampletimes, time)
          } else {
            for (j in 1:samples.per.time) {
              if (length(strain.log[[i]])==1) {
                pickgrp <- strain.log[[i]]
              } else {
                pickgrp <- sample(strain.log[[i]], 1, prob=freq.log[[i]])
              }
              sampleWGS <- c(sampleWGS, pickgrp)
              if (0%in%sampleWGS) {
                stop("Sampled zeroes")
              }
              sampleID <- c(sampleID, i)
              samplepick <- c(samplepick, j)
              sampletimes <- c(sampletimes, time)
            }
          }
          if (samp.schedule != "random" && rec.times[which(inf.ID==i)] > time + samp.freq) {
            sample.times[which(inf.ID==i)] <- time + samp.freq
          }
        }
      }
      # clean up libr etc.
      if (!full) {
        deleters <- NULL
        uniquestrains <- 0
        for (j in 1:length(totcurstrains)) {
          tottype <- 0
          for (k in cur.inf) {
            if (totcurstrains[j]%in%strain.log[[k]]) { # if strain is extant
              tottype <- tottype+1
            }
            if (sum(!strain.log[[k]]%in%totcurstrains)>0) {
              stop("Deleted sequence for observed sample")
            }
          }
          if (tottype>0) { # don't delete if still around
            uniquestrains <- uniquestrains+1
          } else if (tottype==0 && !totcurstrains[j]%in%sampleWGS) { # if not around AND not logged
            deleters <- c(deleters, j) # delete
          }
        }
        if (length(deleters)>0) {
          for (i in sort(deleters, decreasing=T)) {
            libr[[i]] <- NULL
            mut.nuc[[i]] <- NULL
          }
          deletegroup <- totcurstrains[deleters]
          totcurstrains <- totcurstrains[-deleters]
        }
      }
      if (length(cur.sus)==0) {
        eff.cur.inf <- inf.ID[which(sample.times>time)]
      }
      if (length(cur.sus)==0 && sum(sample.times>time)==0) {
        break
      }
    }
    if (tot.inf>mincases) {
      trigger <- TRUE
    } else {
      cat("Insufficient number of infections! (mincases=", mincases, ")\n", sep="")
    }
  }
  if (full) {
    return(invisible(list(epidata=cbind(inf.ID, inf.times, rec.times, inf.source), 
                          sampledata=cbind(pID, sampleID, sampletimes), obs.freq=obs.freq, obs.strain=obs.strain,
                          libr=libr, nuc=mut.nuc, librstrains=totcurstrains, endtime=time)))
  } else {
    return(invisible(list(epidata=cbind(inf.ID, inf.times, rec.times, inf.source), 
                          sampledata=cbind(sampleID, sampletimes, sampleWGS),
                          libr=libr, nuc=mut.nuc, librstrains=totcurstrains, endtime=time)))
  }
}
