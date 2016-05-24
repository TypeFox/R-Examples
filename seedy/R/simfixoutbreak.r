simfixoutbreak <-
  function(ID, inf.times, rec.times, inf.source, mut.rate, equi.pop=10000, shape=flat,
           inoc.size=1, imp.var=25, samples.per.time=1, samp.schedule="random", 
           samp.freq=500, full=FALSE, feedback=500, glen=100000, 
           ref.strain=NULL, ...) {
    
    # WARNINGS
    
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
    if (samples.per.time%%1!=0 || samples.per.time<1) {
      stop("samples.per.time must be a postive integer")
    }
    if (length(inf.times)!=length(rec.times) || length(inf.times)!=length(inf.source)) {
      stop("Infection times, recovery times, infection sources must be vectors of the sample length")
    }
    if (sum(inf.times>=rec.times)>0) {
      stop("Infection times must be earlier than recovery times")
    }
    for (i in 1:length(inf.times)) {
      if (inf.source[i]!=0) {
        if (inf.times[which(ID==inf.source[i])]>=inf.times[i]) {
          stop(paste("Source's infection time after recipient's, p", i, sep=""))
        } else if (rec.times[which(ID==inf.source[i])]<=inf.times[i]) {
          stop(paste("Source's recovery time earlier than recipient's infection time, p", i, sep=""))
        }
      }
    }
    if (mut.rate<0 || mut.rate>=1) {
      stop("Mutation rate must be between 0 and 1")
    }
    if (!is.function(shape)) {
      stop("'shape' must be a function")
    } 
    if (equi.pop%%1!=0 || equi.pop<=0) {
      stop("Equilibrium population size must be a postive integer")
    }
    if (!samp.schedule%in%c("random", "calendar", "individual")) {
      stop("samp.schedule must be 'random', 'calendar', or 'individual'")
    }
    
    #########################
    
    cat("\nSimulating genomic data:\n")
    
    time <- min(inf.times) # in bacterial generations
    newinfect <- 0
    eff.cur.inf <- NULL
    init.inf <- sum(inf.times==min(inf.times))
    init.sus <- length(inf.times)-init.inf
    cur.inf <- ID[which(inf.times==min(inf.times))] # vector of infected person IDs
    n.sus <- length(inf.times)-init.inf
    tot.inf <- init.inf
    sample.times <- numeric(length(inf.times)) # What is the first sampling time for ind. i?
    for (i in 1:length(inf.times)) {
      if (samp.schedule=="random") {
        sample.times[i] <- sample(inf.times[i]:(rec.times[i]-1),1)
      } else if (samp.schedule == "individual") {
        sample.times[i] <- inf.times[i] + samp.freq
      } else if (samp.schedule == "calendar") {
        sample.times[i] <- time+samp.freq
      }
    }
    
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
    w <-0
    for (i in cur.inf) {
      w <- w+1
      nmuts <- rpois(1,imp.var)
      if (nmuts>0) {
        libr[[w]] <- sample(glen,nmuts)
        mut.nuc[[w]] <- numeric(nmuts)
        for (j in 1:nmuts) {
          mut.nuc[[w]][j] <- sample((1:4)[-ref.strain[libr[[w]][j]]], 1)
        }
      }   
      freq.log[[which(ID==i)]] <- 1
      strain.log[[which(ID==i)]] <- w
    }
    
    for (i in (1:length(inf.times))[-which(ID%in%cur.inf)]) {
      freq.log[[i]] <- 0
      strain.log[[i]] <- 0
    }
    
    current.infected <- init.inf
    types <- init.inf # Cumulative number of strain types
    
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
    
    while (time <= max(rec.times)) { # Cycle through bacterial generations until epidemic ceases
      time <- time+1
      if (time%in%rec.times) { # recovery?
        recover <- ID[which(rec.times==time)] # who has recovered?
        cur.inf <- cur.inf[-which(cur.inf%in%recover)] # remove infective(s)
        for (r in 1:length(recover)) {
          strain.log[[which(ID==recover[r])]] <- 0
          freq.log[[which(ID==recover[r])]] <- 0
        }
        if (length(cur.inf)==0) { # If no more infectives
          cat("t=", time, ", S=", n.sus, ", I=", length(cur.inf), 
              ", total genotypes=0\n", sep="")
        }
      }
      if (time%%feedback==0 && length(cur.inf)>0) { # output current status every x generations
        cat("t=", time, ", S=", n.sus, ", I=", length(cur.inf), 
            ", total genotypes=", length(unique(as.numeric(unlist(strain.log)))), ", next rec time=", 
            min(rec.times[which(rec.times>time)]), sep="")
        if (n.sus>0) {
          cat("\n")
        } else {
          cat(", final removal time=", max(rec.times), "\n", sep="")
        }
      }
      
      if (time %in% inf.times) { # infection?
        n.sus <- n.sus-1
        infnow <- ID[which(inf.times==time)]
        for (k in 1:length(infnow)) {
          newinfect <- infnow[k]
          tot.inf <- tot.inf+1
        
          cur.inf <- c(cur.inf, newinfect) # add to current infectives
          ###############
          # pass on strain
          if (inf.source[which(ID==newinfect)]==0) {
            nmuts <- rpois(1,imp.var)
            if (nmuts>0) {
              types <- types+1
              totcurstrains <- c(totcurstrains, types)
              libr[[length(totcurstrains)]] <- sample(glen,nmuts)
              mut.nuc[[length(totcurstrains)]] <- numeric(nmuts)
              for (j in 1:nmuts) {
                mut.nuc[[length(totcurstrains)]][j] <- sample((1:4)[-ref.strain[libr[[length(totcurstrains)]][j]]], 1)
              }
              strain.log[[which(ID==newinfect)]] <- types
              freq.log[[which(ID==newinfect)]] <- 1
            } else {
              strain.log[[which(ID==newinfect)]] <- 1
              freq.log[[which(ID==newinfect)]] <- 1
            }
          } else {
            src <- cur.inf[which(cur.inf==inf.source[which(ID==newinfect)])] # Source of infection
            if (length(strain.log[[which(ID==src)]])==1) { # if source has clonal infection
              inoc.samp <- rep(strain.log[[which(ID==src)]], inoc.size)
              if (0%in%inoc.samp) {
                stop("Zeroes in inoculum")
              }
            } else {
              inoc.samp <- sample(strain.log[[which(ID==src)]], inoc.size, 
                                  prob=freq.log[[which(ID==src)]], replace=T) # take random sample
              if (0%in%inoc.samp) {
                stop("Zeroes in inoculum")
              }
            }
            strain.log[[which(ID==newinfect)]] <- unique(inoc.samp) # distinct types in new infection
            f <- numeric(length(unique(inoc.samp)))
            w <- 1
            for (i in unique(inoc.samp)) {
              f[w] <- sum(inoc.samp==i)
              w <- w+1
            }
            freq.log[[which(ID==newinfect)]] <- f # frequency of types
          }
        }
      }
      # mutate existing strains for each individual
      if (length(cur.inf>0)) {
        #if (is.null(eff.cur.inf)) {
        # cinf <- inf.ID[which(inf.ID%in%cur.inf)]
        #} else {
        #  cinf <- eff.cur.inf
        #}
        for (i in cur.inf) {
          pop.size <- sum(freq.log[[which(ID==i)]])
          #death.prob <- min(0.5 + 0.5*(pop.size-shape(time,span=rec.times[which(ID==i)]-inf.times[which(ID==i)]+1,equi.pop,...))/shape(time,span=rec.times[which(ID==i)]-inf.times[which(ID==i)]+1,equi.pop,...),1)
          death.prob <- min(0.5 + 0.5*(pop.size-shape(time,span=rec.times[which(ID==i)]-inf.times[which(ID==i)]+1,equi.pop))/shape(time,span=rec.times[which(ID==i)]-inf.times[which(ID==i)]+1,equi.pop),1)
          
          if (length(freq.log)==0 || sum(is.na(freq.log))>0 || death.prob>1 || death.prob<0) {
            cat("deathprob=", death.prob, "\npop.size=", pop.size, "\nequi.pop=", equi.pop, "\nFreq.log:\n")
            if (length(freq.log)>0) {
              for (k in 1:length(freq.log)) {
                cat(freq.log[[which(ID==i)]][k], "\n")
              }
            }
          }
          freq.log[[which(ID==i)]] <- 2*rbinom(length(freq.log[[which(ID==i)]]), freq.log[[which(ID==i)]], 1-death.prob)
          if (0 %in% freq.log[[which(ID==i)]]) {
            zeros <- which(freq.log[[which(ID==i)]]==0)
            if (length(zeros)==length(freq.log[[which(ID==i)]])) {
              freq.log[[which(ID==i)]] <- 1
              strain.log[[which(ID==i)]] <- strain.log[[which(ID==i)]][1]
            } else {
              freq.log[[which(ID==i)]] <- freq.log[[which(ID==i)]][-zeros]
              strain.log[[which(ID==i)]] <- strain.log[[which(ID==i)]][-zeros]
            }
          }
          if (length(strain.log[[which(ID==i)]])!=length(freq.log[[which(ID==i)]])) {
            stop("Error")
          }
          n.mutations <- rbinom(1, sum(freq.log[[which(ID==i)]]), mut.rate)
          if (n.mutations > 0) {
            for (mt in 1:n.mutations) {
              types <- types+1
              if (length(strain.log[[which(ID==i)]])==1) {
                mutate.grp <- strain.log[[which(ID==i)]]
              } else {
                mutate.grp <- sample(strain.log[[which(ID==i)]], 1, prob=freq.log[[which(ID==i)]])
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
                strain.log[[which(ID==i)]] <- c(strain.log[[which(ID==i)]], types)
                freq.log[[which(ID==i)]] <- c(freq.log[[which(ID==i)]], 1)
                totcurstrains <- c(totcurstrains, types)
              }
            }
          }
        }
      }
      # take samples, make observations
      if (time%in%sample.times) {
        smpat <- which(sample.times==time & rec.times > time & inf.times <= time)
        for (i in smpat) {
          if (full) {
            n <- length(obs.freq)+1
            obs.freq[[n]] <- freq.log[[i]]
            obs.strain[[n]] <- strain.log[[i]]
            sampleID <- c(sampleID, n)
            pID <- c(pID, ID[i])
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
              sampleID <- c(sampleID, ID[i])
              samplepick <- c(samplepick, j)
              sampletimes <- c(sampletimes, time)
            }
          }
          if (sample.times[i]+samp.freq<rec.times[i] && samp.schedule!="random") {
            sample.times[i] <- sample.times[i]+samp.freq
          }
        }
      }
      # clean up libr etc.
      if (!full && length(cur.inf)>0) {
        deleters <- NULL
        uniquestrains <- 0
        for (j in 1:length(totcurstrains)) {
          tottype <- 0
          for (k in cur.inf) {
            if (totcurstrains[j]%in%strain.log[[which(ID==k)]]) { # if strain is extant
              tottype <- tottype+1
            }
            if (sum(!strain.log[[which(ID==k)]]%in%totcurstrains)>0) {
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
          for (i in sort(deleters, decreasing=TRUE)) {
            libr[[i]] <- NULL
            mut.nuc[[i]] <- NULL
          }
          deletegroup <- totcurstrains[deleters]
          totcurstrains <- totcurstrains[-deleters]
        }
      }
      #if (length(cur.sus)==0) {
      #  eff.cur.inf <- which(sample.times>time)
      #}
    }

    if (full) {
      return(invisible(list(epidata=cbind(ID, inf.times, rec.times, inf.source), 
                            sampledata=cbind(pID, sampleID, sampletimes), obs.freq=obs.freq, obs.strain=obs.strain,
                            libr=libr, nuc=mut.nuc, librstrains=totcurstrains, endtime=time)))
    } else {
      return(invisible(list(epidata=cbind(ID, inf.times, rec.times, inf.source), 
                            sampledata=cbind(sampleID, sampletimes, sampleWGS),
                            libr=libr, nuc=mut.nuc, librstrains=totcurstrains, endtime=time)))
    }
  }

