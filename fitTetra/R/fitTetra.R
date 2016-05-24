fitTetra <-
function (
            marker, data, diplo=NA,
            select=TRUE, diploselect=TRUE, # by default all samples of tetraploids and diploids selected
            maxiter=40, maxn.bin=200, nbin=200,
            sd.threshold=0.1, p.threshold=0.99,
            call.threshold=0.6, peak.threshold=0.85,
            try.HW=TRUE, dip.filter=1,
            sd.target=NA,
            plot="none", plot.type="png", plot.dir=NA) {
  if (class(data)!="data.frame" || 
    length(which(names(data)=="MarkerName"))!=1 ||
    length(which(names(data)=="SampleName"))!=1 ||
    length(which(names(data)=="ratio"))!=1) {
    stop("data is not a valid data frame")
  }  
  # dip.filter: change to integer if needed 
  # (for compatibility with earlier versions where it was logical)          
  dip.filter <- as.integer(dip.filter)
  if (dip.filter<0 || dip.filter>2) dip.filter <- 1
  # arrange plotting output:
  plot <- tolower(plot)
  if (plot!="none") {
    if (substr(plot.type,1,1)=="_") 
      plot.type <- substr(plot.type,2,5) #already checked by saveMarkerModels
    else plot.type <- check.plottype(plot.type)
    if (plot.type=="none") plot<-"none"
  }
  plotall <- plot=="all"
  plotfitted <- plotall || plot=="fitted"
  if (is.na(plot.dir)) plot.dir <- ""
  else {
    plot.dir <- gsub("(^ +)|( +$)", "", plot.dir) #strip leading and trailing blanks
    if (plot.dir != "") {
      plot.dir <- paste(plot.dir,"/",sep="") 
    }
  }
  # get all markernames, samplenames and diplonames:
  markernames <- levels(as.factor(as.character(data$MarkerName)))
  markername <- markernames[marker]
  mrknr <- padded(marker,length(markernames)) #adds leading zeroes if needed
  samplenames <- levels(as.factor(as.character(data$SampleName)))
  if (class(diplo)=="data.frame") #not NA
    diplonames <- levels(as.factor(as.character(diplo$SampleName)))
  # start log
  log <- ""
  log <- append(log,paste(marker,"\tname=",markername,sep=""))
  resultnames <- getResultnames(ng=5)
  rejected <- FALSE
  
  # first: select the diploid data (if any) to plot into the tetraploid histogram:
  drr <- numeric(0) #for later test on length>0
  if (class(diplo)=="data.frame") {
    if (class(diplo$ratio)!="numeric") stop("diplo$ratio is not numeric")
    if (length(diploselect)!=nrow(diplo)) {
      diploselect <- rep(diploselect, times=ceiling(nrow(diplo)/length(diploselect)))
      if (length(diploselect)>nrow(diplo)) diploselect <- diploselect[1:nrow(diplo)]
    }
    drr <- diplo$ratio[(diplo$MarkerName==markername) & diploselect & !is.na(diplo$ratio)]
    #drr contains only non-NA ratios and excludes any for which diploselect is FALSE
  }

  #now the tetraploids:
  if (class(data$ratio)!="numeric") stop("data$ratio is not numeric")
  if (length(select)!=nrow(data)) {
    select <- rep(select, times=ceiling(nrow(data)/length(select)))
    if (length(select)>nrow(data)) select <- select[1:nrow(data)]
  }
  nsamp<-length(samplenames)
  rr <- data$ratio[(data$MarkerName==markername) & select & !is.na(data$ratio)] #selection for minimum intensity
  #rr contains only non-NA ratios and excludes any for which select is FALSE
  # calculate the total number of selected samples, including the ones with ratio NA even if they are not in data:
  missingsamp <- nsamp - length(data$ratio[data$MarkerName==markername])
  countselected <- sum(select[data$MarkerName==markername], na.rm=T)
  nselsamp <- countselected + missingsamp; rm(missingsamp,countselected)
  if (length(rr)<50) { # in CodomMarker at least 10*ng observations required
    rejected <- TRUE
    log <- append(log,paste(marker,"","","discarded: too few data",sep="\t"))
    mrkresult <- NA
  }
  else {
    if (length(rr)<nselsamp*call.threshold) {
      rejected <- TRUE
      log <- append(log,paste(marker,"","","discarded: too many NA ratios",sep="\t"))
      mrkresult <- NA
    }
    else { #test if no samples occur more than once for this marker
      allsamp <- as.factor(data$SampleName[data$MarkerName==markername])
      tabsamp <- tabulate(allsamp)
      if (max(tabsamp)>1) {
        rejected <- TRUE
        log <- append(log,paste(marker,"","","discarded: some samples occur more than once",sep="\t"))
        mrkresult <- NA
      }
      if (class(diplo)=="data.frame") {
        allsamp <- diplo$SampleName[diplo$MarkerName==markername]
        tabsamp <- tabulate(allsamp)
        if (max(tabsamp)>1) {
          rejected <- TRUE
          log <- append(log,paste(marker,"","","discarded: some diploid samples occur more than once",sep="\t"))
          mrkresult <- NA
        }
      }  
    }
  }  

  if (!rejected) { 
    # sufficient observations and no double samples 
    names(rr) <- data$SampleName[(data$MarkerName==markername) & select & !is.na(data$ratio)] #tetraploid sample names, should be unique within marker
    if (class(diplo)=="data.frame") {
      names(drr) <- diplo$SampleName[(diplo$MarkerName==markername) & diploselect & !is.na(diplo$ratio)] #diploid sample names, should be unique within marker
    }
    # try several models with 5 peaks
    ng <- 5 #ng = number of genotypic classes
    nmodel <- 32 #8 models with 2, 3 or 4 start configurations
    mrkresult <- list()
    mrkresult$stats<-data.frame(rep(NA,nmodel))
    for (i in 2:length(resultnames)) mrkresult$stats[[i]] <- rep(NA,nmodel)
    names(mrkresult$stats) <- resultnames
    mrkresult$stats$m <- 1:nmodel
    mrkresult$probs <- list()

    # 8 models with clustering:
    plot.fname <- paste(plot.dir,"plots ",mrknr," A ",markername,sep="")
    if (plotall) start.plot(plot.type,plot.fname, ncol=2, nrow=4) 
    # use kmeans with fixed random seed for reproducible results
    # but don't affect the random or non-random number sequence for
    # the calling program:
    if (!exists('.Random.seed')) runif(1) #.Random.seed only generated in first call of random number generator
    savedseed <- .Random.seed; set.seed(3)

    # initialise the cluster means once
    clusinit <- ClusterInit(asin(sqrt(rr)), ng=ng)   # returns clus.mu and clus.sd on transformed scale
    clus.mu <- sin(clusinit$clus.mu)^2  # transform clus.mu back to 0-1 scale, but keep clus.sd on transformed scale
    
    for (model in 1:8) if (try.HW || getPtype(model)!="p.HW") {
      mrkresult<-MarkerResult(marker,markername,rr,model=model,mutype=getMutype(model),ptype=getPtype(model),
                              maxiter=maxiter,maxn.bin=maxn.bin,nbin=nbin, mustart=clus.mu, sdstart=clusinit$clus.sd,
                              mrkresult=mrkresult,nsamp=nselsamp,plothist=plotall)
      if (plotall && is.na(mrkresult$stats$npar[model]))
        errorplot(paste(marker,markername,getModelName(getMutype(model),getPtype(model))))
    }
    .Random.seed <- savedseed
    if (plotall) save.plot(plot.type,plot.fname)
    #same models without clustering:
    plot.fname <- paste(plot.dir,"plots ",mrknr," B ",markername,sep="")
    if (plotall) start.plot(plot.type,plot.fname, ncol=2, nrow=4)
    for (model in 9:16) if (try.HW || getPtype(model)!="p.HW") {
      mrkresult<-MarkerResult(marker,markername,rr,model=model,mutype=getMutype(model),ptype=getPtype(model),clus=F,
                              maxiter=maxiter,maxn.bin=maxn.bin,nbin=nbin,mrkresult=mrkresult,nsamp=nselsamp,plothist=plotall)
      if (plotall && is.na(mrkresult$stats$npar[model]))
        errorplot(paste(marker,markername,getModelName(getMutype(model),getPtype(model))))
    }
    if (plotall) save.plot(plot.type,plot.fname)
    modelsfitted <- 16
    mrkresult$stats <- calcSelcrit(mrkresult$stats,rep(T,nrow(mrkresult$stats)),sd.target)
    optmodel1 <- which.min(mrkresult$stats$selcrit)
    if (length(optmodel1)==0 || is.na(optmodel1)) {
      rejected <- TRUE
      log <- append(log,paste(marker,"optmodel1=NA","","discarded: optmodel1=NA",sep="\t"))
    } else {
      #optmodel1!=NA
      #First we check for two situations indication a possibly wrong fit:
      # - the nulli or quadruplex peak is very small while the simplex or triplex is big, or
      # - there is a small peak at the simplex or triplex while the flanking peaks are bigger
      #In these cases we run all the models again with a new mustart
      opt.ptype <- getPtype(optmodel1)
      opt.mutype <- getMutype(optmodel1)
      logline <- paste(marker,"\toptmodel1=",optmodel1,"\t",getModelName(opt.mutype,opt.ptype),sep="")
      #Check for a small nulli or quadruplex next to a large simplex or triplex:
      #Note: in F1 sim+duplex or du+triplex are also possible, but nulli+sim
      #and tri+quadri are frequently shifted towards the middle
      P0 <- which(names(mrkresult$stats)=="Pact0")
      Pactpeak <- unlist(subset(mrkresult$stats, subset=1:nrow(mrkresult$stats)==optmodel1)[P0:(P0+4)])
      #Pactpeak[1:5] are now the fractions of selected samples in nulliplex-tetraplex peaks
      max.nulli <- 0.025; min.sim<-0.15 #thresholds for small extreme peak detection
                                        #false positive in F1 at 1:8:18:8:1 (1 of 25 seg.types)
                                        #and in HW with major allele freq 0.58-0.61
      if ( (Pactpeak[1]<max.nulli && Pactpeak[2]>min.sim) ||
           (Pactpeak[5]<max.nulli && Pactpeak[4]>min.sim) ) {
        loside <- Pactpeak[1]<max.nulli && Pactpeak[2]>min.sim #at low side of range
        if (loside) {
          if (Pactpeak[5]<max.nulli && Pactpeak[4]>min.sim) {
            #situation occurs at both sides, which is most extreme?
            loside <- (Pactpeak[2]-Pactpeak[1]) > (Pactpeak[4]-Pactpeak[5])
          }
        }
       
        #set the nwmu according to low side or high side:  
        nwmu <- numeric(5)
        mu0 <- which(names(mrkresult$stats)=="muact0")
        muactpeak <- unlist(subset(mrkresult$stats, subset=1:nrow(mrkresult$stats)==optmodel1)[mu0:(mu0+4)]) #the current mu's
        if (loside) {
          nwmu[1:4] <- muactpeak[2:5]
          nwmu[5] <- (nwmu[4]+pi/2)/2
        } else {
          nwmu[2:5] <- muactpeak[1:4]
          nwmu[1] <- nwmu[2]/2
        }  
        #fit all models with the new mustart:
        nwmu <- sin(nwmu)*sin(nwmu) #to untransformed ratios
        plot.fname <- paste(plot.dir,"plots ",mrknr," C ",markername,sep="")
        if (plotall) start.plot(plot.type,plot.fname, ncol=2, nrow=4)
        for (model in 17:24) if (try.HW || getPtype(model)!="p.HW") {
          mrkresult<-MarkerResult(marker,markername,rr,model=model,mutype=getMutype(model),ptype=getPtype(model),clus=F,
                                  maxiter=maxiter,maxn.bin=maxn.bin,nbin=nbin,mustart=nwmu,
                                  mrkresult=mrkresult,nsamp=nselsamp,plothist=plotall)
          if (plotall && is.na(mrkresult$stats$npar[model]))
            errorplot(paste(marker,markername,getModelName(getMutype(model),getPtype(model))))
        }
        if (plotall) save.plot(plot.type,plot.fname)
        modelsfitted <- 24 
      } # small nulli or quadruplex peak next to big peak
      else {
        #check for dip at simplex or triplex peak:
        if (mrkresult$stats$dip[optmodel1]) {
          # (note: indices from [1] to [5], not [0] to [4] !)
          p<-rep(0,5)
          p[1]<-mrkresult$stats$P0[optmodel1]
          p[2]<-mrkresult$stats$P1[optmodel1]
          p[3]<-mrkresult$stats$P2[optmodel1]
          p[4]<-mrkresult$stats$P3[optmodel1]
          p[5]<-mrkresult$stats$P4[optmodel1]
          mubk<-rep(0,5)
          mubk[1]<-mrkresult$stats$mu0[optmodel1]
          mubk[2]<-mrkresult$stats$mu1[optmodel1]
          mubk[3]<-mrkresult$stats$mu2[optmodel1]
          mubk[4]<-mrkresult$stats$mu3[optmodel1]
          mubk[5]<-mrkresult$stats$mu4[optmodel1]
          gsamp<-hist(mrkresult$probs[[optmodel1]]$maxgeno,breaks=c(-0.5,0.5,1.5,2.5,3.5,4.5),plot=F)$counts #based on maxgeno: all samples in peaks
          toler<-1.02 #doesn't harm to do this too often, we still should obtain the correct model
          #test both valleys and select widest (we do not test for valley at duplex peak)
          valley<-numeric(0)
          if ( (gsamp[2]<(toler*gsamp[1]) && gsamp[2]<(toler*max(gsamp[3:5])) && gsamp[1]>3)  ||
               (opt.ptype=="p.free" && p[2]<(toler*p[1]) && p[2]<(toler*max(p[3:5])) && p[1]>0.01) ) {
             valley <- 1
          }
          if ( (gsamp[4]<(toler*gsamp[5]) && gsamp[4]<(toler*max(gsamp[1:3])) && gsamp[5]>3 ) ||
                    (opt.ptype=="p.free" && p[4]<(toler*p[5]) && p[4]<(toler*max(p[1:3])) && p[5]>0.01) ) {
             valley[length(valley)+1] <- 3
          }
          if (length(valley)==2) {
            #select widest valley based on the mu's on the transformed scale
            if ( (mrkresult$stats$mu4[optmodel1]-mrkresult$stats$mu2[optmodel1]) >
                 (mrkresult$stats$mu2[optmodel1]-mrkresult$stats$mu0[optmodel1]) ) {
              valley <- 3
            }
            else valley <- 1
            logline <- paste(logline,"\tdips at peak 1&3, widest at",valley,sep="")
          }
          else if (length(valley)==1)
            logline <- paste(logline,"\tdip at peak",valley,sep="")
          if (length(valley)==1) {
            #calculate the new mustart:
            if (valley==1) {
              mubk[2:4]<-mubk[3:5]
              mubk[5]<-(mubk[5]+1)/2
            }
            else { #valley==3
              mubk[2:4]<-mubk[1:3]
              mubk[1]<-mubk[1]/2
            }
            # fit all models with the new mustart:
            plot.fname <- paste(plot.dir,"plots ",mrknr," C ",markername,sep="")
            if (plotall) start.plot(plot.type,plot.fname, ncol=2, nrow=4)
            for (model in 17:24) if (try.HW || getPtype(model)!="p.HW") {
              mrkresult<-MarkerResult(marker,markername,rr,model=model,mutype=getMutype(model),ptype=getPtype(model),clus=F,
                                      maxiter=maxiter,maxn.bin=maxn.bin,nbin=nbin,mustart=mubk,mrkresult=mrkresult,nsamp=nselsamp,plothist=plotall)
              if (plotall && is.na(mrkresult$stats$npar[model]))
                errorplot(paste(marker,markername,getModelName(getMutype(model),getPtype(model))))
            }
            if (plotall) save.plot(plot.type,plot.fname)
            modelsfitted <- 24
            #check if the best model still has a dip:
            mrkresult$stats <- calcSelcrit(mrkresult$stats,rep(T,nrow(mrkresult$stats)),sd.target) #all rows
            optmodel2 <- which.min(mrkresult$stats$selcrit)
            if (mrkresult$stats$dip[optmodel2]) {
              #optimal model still has dip 
              #we again take the mu's and the valley of optmodel1 but now shift in opposite direction:
              mubk[1]<-mrkresult$stats$mu0[optmodel1]
              mubk[2]<-mrkresult$stats$mu1[optmodel1]
              mubk[3]<-mrkresult$stats$mu2[optmodel1]
              mubk[4]<-mrkresult$stats$mu3[optmodel1]
              mubk[5]<-mrkresult$stats$mu4[optmodel1]
              #calculate the new mustart:
              if (valley==1) {
                mubk[2]<-mubk[1]
                mubk[1]<-mubk[1]/2
              }
              else { #valley==3
                mubk[4]<-mubk[5]
                mubk[5]<-(mubk[5]+1)/2
              }
              # fit all models with the new mustart:
              plot.fname <- paste(plot.dir,"plots ",mrknr," D ",markername,sep="")
              if (plotall) start.plot(plot.type,plot.fname, ncol=2, nrow=4)
              for (model in 25:32) if (try.HW || getPtype(model)!="p.HW") {
                mrkresult<-MarkerResult(marker,markername,rr,model=model,mutype=getMutype(model),ptype=getPtype(model),clus=F,
                                        maxiter=maxiter,maxn.bin=maxn.bin,nbin=nbin,mustart=mubk,mrkresult=mrkresult,nsamp=nselsamp,plothist=plotall)
                if (plotall && is.na(mrkresult$stats$npar[model]))
                  errorplot(paste(marker,markername,getModelName(getMutype(model),getPtype(model))))
              }
              if (plotall) save.plot(plot.type,plot.fname)
              modelsfitted <- 32
            } # if dip[optmodel2]
          } #if length(valley)==1
        } # if dip[optmodel1]
      } # 2 big peaks else 
      log <- append(log, logline)
      #now select the best model (optmodel2) taking dip into account:
      if (dip.filter>0) {
        dipok<-!(is.na(mrkresult$stats$dip) | mrkresult$stats$dip) #only rows with model without dip
        mrkresult$stats <- calcSelcrit(mrkresult$stats,dipok,sd.target)
        optmodel2 <- which.min(mrkresult$stats$selcrit[1:modelsfitted])
        if (length(optmodel2)==0 || is.na(optmodel2)) {
          mrkresult$stats <- calcSelcrit(mrkresult$stats,rep(T,nrow(mrkresult$stats)),sd.target) #all rows
          optmodel2 <- which.min(mrkresult$stats$selcrit[1:modelsfitted])
        }
      }
      else {
        mrkresult$stats <- calcSelcrit(mrkresult$stats,rep(T,nrow(mrkresult$stats)),sd.target) #all rows
        optmodel2 <- which.min(mrkresult$stats$selcrit[1:modelsfitted])
      } 
      if (length(optmodel2)==0 || is.na(optmodel2)) {
        # discard this marker
        rejected <- TRUE
        log <- append(log,paste(marker,"optmodel2=NA","","marker discarded: optmodel2=NA",sep="\t"))
      } else {
        # optmodel2 found
        opt.ptype <- getPtype(optmodel2)
        opt.mutype <- getMutype(optmodel2)
        logline <- paste(marker,"\toptmodel2=",optmodel2,"\t",getModelName(opt.mutype,opt.ptype),sep="")
        if (mrkresult$stats$sdtrans0[optmodel2]>sd.threshold) {
          rejected <- TRUE
          log <- append(log,paste(logline,"discarded: sd>sd.threshold",sep="\t"))
          mrkresult$stats$message[optmodel2] <- "rejected: sd>sd.threshold"
        } 
        else if (mrkresult$stats$dip[optmodel2] && dip.filter==2) {
          rejected <- TRUE
          mrkresult$stats$message[optmodel2] <- "rejected: no models without dip"
          log <- append(log,paste(logline, "discarded: no models without dip"))
        }
        else {
          # geno<-maxgeno for samples above p.threshold, else NA:
          mrkresult$probs[[optmodel2]]$geno<-rep(NA,length(rr))
          sel.maxP <- mrkresult$probs[[optmodel2]]$maxP
          sel.maxP <- sel.maxP>=p.threshold                       #without the temporary fix, not needed any more                 
          mrkresult$probs[[optmodel2]]$geno[sel.maxP] <- mrkresult$probs[[optmodel2]]$maxgeno[sel.maxP]
          if (length(mrkresult$probs[[optmodel2]]$geno[!is.na(mrkresult$probs[[optmodel2]]$geno)])<call.threshold*nselsamp) {
              # note that we compare to the samples with select=T, not to the total number of samples
              rejected <- TRUE
              log <- append(log,paste(logline,"discarded: less than minimum number of samples called",sep="\t"))
              mrkresult$stats$message[optmodel2] <- "rejected: less than minimum number of samples called"
          } else {
            gsamp<-hist(mrkresult$probs[[optmodel2]]$geno,breaks=c(-0.5,0.5,1.5,2.5,3.5,4.5),plot=F)$counts   #based on called samples
            if (max(gsamp)/sum(gsamp)>peak.threshold) {
              rejected <- TRUE
              log <- append(log,paste(logline, "discarded: more than maximum fraction of samples in one peak",sep="\t"))
              mrkresult$stats$message[optmodel2] <- "rejected: more than maximum fraction of samples in one peak"
            } else {
              #optmodel2 accepted
              if (mrkresult$stats$dip[optmodel2] && dip.filter==1) {
                #if dip.filter==2 marker was already rejected; if 0 there might be a model without dip but worse selcrit
                logline <- paste(logline, "no models without dip")
              }  
              log <- append(log,logline)
              probs <- makeprobs(marker,markername,samplenames,getModelName(opt.mutype,opt.ptype),select,
                                 data$MarkerName==markername,rr,mrkresult$probs[[optmodel2]],data)
              if (plotfitted) {
                #plot the fitted model with the assigned genotypes:
                modeldata <- list()
                modeldata$psi <- list()
                nfirst <- which(names(mrkresult$stats)=="mutrans0")
                modeldata$psi$mu <- as.numeric(mrkresult$stats[optmodel2,nfirst:(nfirst+ng-1)])
                modeldata$psi$sigma <- as.numeric(mrkresult$stats[optmodel2,(nfirst+ng):(nfirst+2*ng-1)])
                modeldata$psi$p <- as.numeric(mrkresult$stats[optmodel2,(nfirst+2*ng):(nfirst+3*ng-1)])
                plotfitted.fname <- paste(plot.dir,mrknr," ",markernames[marker],sep="") #final plot, fitted
                start.plot(plot.type,plotfitted.fname)
                layout(matrix(c(1,2)), heights=c(2,1)) #2 plots, upper 3x as high as lower
                #upper plot: the fitted model
                par(mar=c(0,4,3,1))
                htet <- PlotHistDensity(rr, modeldata, xlim=c(0.001,0.999), trafo="asr",
                                                maintitle=paste(marker,markername),
                                                subtitle="", xlabel="", xaxis="n",
                                                nbreaks=40, frequ=T)
                if (length(drr)>0) {
                  #plot diploids histogram into tetraploids plot:
                  hdip <- hist(drr,breaks=htet$breaks,plot=F)
                  hdip$counts[hdip$counts>htet$ylim[2]]<-htet$ylim[2]
                  nbreaks <- length(htet$breaks)-1
                  binwidth <- 1/nbreaks
                  barwidth <- 1/3 #relative to binwidth
                  par(new=T)
                  barplot (hdip$counts,
                    width=binwidth*barwidth,
                    space=1/barwidth-1,
                    xlim=c(binwidth*(1-barwidth)/2,1+binwidth*(1-barwidth)/2),
                    ylim=htet$ylim,col="gray",main="",sub="",xlab="",xaxt="n")
                }
                #lower plot: the assigned genotypes
                mrkresult$probs[[optmodel2]]$geno[is.na(mrkresult$probs[[optmodel2]]$geno)] <- -2 #plot well separated from the valid scores 0:4
                par(mar=c(5,4,0,1))
                sampcol = rep("blue",length(rr))
                sampcol[mrkresult$probs[[optmodel2]]$geno < 0] = "red"
                plot(0:1~0:1, type="n",
                  ylab="geno",ylim=c(-3,5),xlab="signal ratio",xlim=c(0.001,0.999) )
                lines(c(-2,-2)~c(par("usr")[1]+0.004,par("usr")[2]-0.004),col="lightgrey",lty=3)
                lines(c(0,0) ~ c(par("usr")[1]+0.004,par("usr")[2]-0.004),col="lightgrey",lty=3)
                lines(c(1,1) ~ c(par("usr")[1]+0.004,par("usr")[2]-0.004),col="lightgrey",lty=3)
                lines(c(2,2) ~ c(par("usr")[1]+0.004,par("usr")[2]-0.004),col="lightgrey",lty=3)
                lines(c(3,3) ~ c(par("usr")[1]+0.004,par("usr")[2]-0.004),col="lightgrey",lty=3)
                lines(c(4,4) ~ c(par("usr")[1]+0.004,par("usr")[2]-0.004),col="lightgrey",lty=3)
                points(mrkresult$probs[[optmodel2]]$geno~rr,
                  pch="|",cex=0.5,col=sampcol)
                save.plot(plot.type,plotfitted.fname)
              } # if plotfitted
            } # optmodel2 accepted (not too many samples in one peak)
          }  # not too few sample genotypes called
        } # optmodel2 not rejected for large sigma / BIC>0
      } #optmodel2!=NA
    } #optmodel1!=NA
  } # if sufficient observations
  if (rejected) {
    if (plotfitted) {
      plotfitted.fname <- paste(plot.dir,"rejected_",mrknr," ",markernames[marker],sep="") #final plot, unfitted
      start.plot(plot.type,plotfitted.fname)
      layout(matrix(c(1,1)))
      ymax <- 40
      nbreaks <- 40
      par(mar=c(5,4,3,1))
      if (length(rr[!is.na(rr)])>0) {
        nbars <- min(ceiling(nbreaks/(max(rr)-min(rr))),100,ceiling(length(rr)/2))
        if (nbars<nbreaks) nbars<-nbreaks
        htet <- hist(rr,breaks=seq(0,1,by=1/nbars),plot=FALSE)
        if (ymax<1.1*max(htet$counts)) ymax<-1.1*max(htet$counts)
      } else {
        # hist fails without data, so we create a custom histogram:
        htet <- list(breaks=seq(0,1,by=1/nbreaks),counts=rep(0,nbreaks))
        class(htet) <- "histogram"
      }
      #plot the tetraploids histogram:
      barplot (htet$counts,
            width=htet$breaks[2]-htet$breaks[1],
            space=0,xlim=c(0.001,0.999),ylim=c(0,ymax), main=paste(marker,markername,"(no fit)"),
            col="white", xlab="signal ratio",ylab="Frequency")
      axis(1,labels=TRUE, at=seq(0,1, by=0.2))
      if (length(drr)>0) {
        #plot diploids histogram into tetraploids plot:
        hdip <- hist(drr,breaks=htet$breaks,plot=F)
        hdip$counts[hdip$counts>ymax]<-ymax
        binwidth <- 1/length(hdip$counts)
        barwidth <- 1/3 #relative to binwidth
        par(new=T)
        barplot (hdip$counts,                                                
          width=binwidth*barwidth,
          space=1/barwidth-1,
          xlim=c(binwidth*(1-barwidth)/2,1+binwidth*(1-barwidth)/2),
          ylim=c(0,ymax),
          col="gray",main="",sub="",xlab="",xaxt="n")
      }
      save.plot(plot.type,paste("rejected_",plotfitted.fname,sep=""))
    } #plotfitted
  } #rejected
  # finally prepare the output list:
  result <- list(log=log)
  if (length(mrkresult)==1) { #only if no fitting attempted because too few samples
    result$allmodeldata <- NA
    result$scores <- NA
    result$diploscores <- NA
    result$modeldata <- data.frame(
      marker=marker,
      markername=markername,
      m=NA,
      model="none",
      nsamp=nselsamp,
      nsel=length(rr))
    for (i in 7:length(resultnames)) result$modeldata[[i]] <- NA
    names(result$modeldata) <- resultnames
    result$modeldata$message <- "rejected: too few observations selected"
  }  
  else {
    result$allmodeldata <- subset(mrkresult$stats,subset=c(rep(T,modelsfitted),rep(F,32-modelsfitted)))
    # Published version: modeldata has only some columns from allmodeldata.
    # New version 6-MAR-2012: modeldata has all columns from allmodeldata
    if (rejected) {
      if (!exists("optmodel2") || length(optmodel2)!=1 || is.na(optmodel2)) {
        result$modeldata <- subset(mrkresult$stats,subset= (1:nrow(mrkresult$stats))==1)
        result$modeldata$m <- NA  #no model selected
        result$modeldata$model <- NA  #no model selected
        if (!exists("optmodel2") || length(optmodel2)!=1 || is.na(optmodel2)) {
          result$modeldata$message <- "optmodel1==NA"
        } else {
          result$modeldata$message <- "optmodel2==NA"
        }  
      } else {
        result$modeldata <- subset(mrkresult$stats,subset= (1:nrow(mrkresult$stats))==optmodel2)
      }    
      for (i in 7:(length(resultnames)-1)) result$modeldata[[i]] <- NA #leave columns up to nsel and message
      result$scores <- NA
      result$diploscores <- NA
    }  
    else { #not rejected
      result$modeldata <- subset(mrkresult$stats,subset= (1:nrow(mrkresult$stats))==optmodel2)
      result$scores <- probs
      if (class(diplo)=="data.frame") { #not NA
        # if diplo exists, also export scores for all samples in diplo:
        psi <- list(mu=c(result$modeldata$mutrans0,result$modeldata$mutrans1,result$modeldata$mutrans2,result$modeldata$mutrans3,result$modeldata$mutrans4),
                sigma=c(result$modeldata$sdtrans0,result$modeldata$sdtrans1,result$modeldata$sdtrans2,result$modeldata$sdtrans3,result$modeldata$sdtrans4),
                p=c(result$modeldata$P0,result$modeldata$P1,result$modeldata$P2,result$modeldata$P3,result$modeldata$P4))
        diploprobs <- EMGaussExp.vectorized(asin(sqrt(drr)), psi) #matrix, per diplo sample ng numbers: probs that sample belongs to each peak
        #calculate maxgeno, maxP, geno:  
        maxPna <- apply(diploprobs,1,max) #find the maxP of each row of p-values (= of each sample)
        maxP <- maxPna[!is.na(maxPna)] # omit missing values
        probs <- data.frame(diploprobs)
        rownames(probs) <- names(drr) #note: does not include samples rejected because select=F
        names(probs) <- c("P0","P1","P2","P3","P4")
        #list for each sample the maximum p value and the maxgeno corresponding to it:
        probs$maxgeno<-max.col(diploprobs)-1  # for every row (sample) the max peak (0..4 instead of 1..5)
        probs$maxP<-maxPna # for every row the P at the max
        #add extra columns and add rows for samples not in drr: 
        result$diploscores <- makeprobs(marker,markername,diplonames,getModelName(opt.mutype,opt.ptype),diploselect,diplo$MarkerName==markername,drr,probs,diplo)
        result$diploscores$geno<-rep(NA,nrow(result$diploscores))
        sel.maxP <- which(result$diploscores$maxP>=p.threshold)
        result$diploscores$geno[sel.maxP] <- result$diploscores$maxgeno[sel.maxP]
      }
      else result$diploscores <- NA
    }  
  }  
  result
}
