plot.IBD<-function(x,kinshipth=NULL, ellipse.coverage=.95,...) {
  # x: an IBD object returned by IBDcheck
  # kinshipth: for plots intended to assign relationships, plot only
  #            pairs with estimated kinship > kinshipth
  # optional arguments passed to plot()

  #vectors of ids
  nsub=nrow(x$snp.data)
  subids<-x$subject.support$subids
  idsub<-pairIDs(nsub,subids) 
  out<-x$ibd.study
  datall <- data.frame(member1=idsub[,1],member2=idsub[,2],
                       pz0=out$pz0,pz1=out$pz1)
  num.allpairs<-nrow(datall)
  
  simulate<-x$simparams$simulate

  # Create an indicator to subset plot of all subjects if kinshipth is non-null
  if(!is.null(kinshipth)) {
    kicf=out$pz2/2+out$pz1/4 #kinship coefficient estimate
    if(kinshipth=="empirical") { # Use 99th %ile of kinship in unrelated
      if(is.null(x$ibd.ur)) {
        stop("Matrix of IBD coefficients for unrelated pairs required when setting kinshipth=`empirical'")
      }
      kiur=x$ibd.ur$pz2/2+x$ibd.ur$pz1/4 
      kinshipth<-quantile(kiur,probs=.99)
    }
    k.ind<-(kicf>kinshipth)
    datall<-datall[k.ind,]
  } 
  
  if(simulate==FALSE){  
    pp <- flagAllPairs(datall,simulate,...) 
  } else {
    if(is.null(x$ibd.ur) && is.null(x$ibd.po) && is.null(x$ibd.fs) && is.null(x$ibd.hs) && is.null(x$ibd.co) && is.null(x$ibd.user)) {
      stop(paste("Couldn't find estimated IBD proportions for any simulated\n",
                 "relationships. Run IBDcheck with `simulate=TRUE`"))
    }

    dat<-data.frame(member1=idsub[,1],member2=idsub[,2],pz0=out$pz0,pz1=out$pz1)
     
    #plot for all pairs includes ellipse based on simulated unrelated pairs,
    #whose level needs a Bonferroni-like adjustment
    coverage.u<-coverage.unrel(ellipse.coverage,num.allpairs)
    flagAllPairs(datall,simulate,x$ibd.ur,coverage.u,...)
    # need to pause between plots, or we won't see this one
    pausePlotting() #kludge

    rships<-c("mz","po","fs","hs","co","user")
    rships.full<-c("MZtwins/duplicates","parent-offspring","full sibs","half sibs","cousins","user")
    rcols<-c("purple","blue","cyan","green","red","yellow")
    pp<-list()
    for(rr in 1:length(rships)) {
      ibd.obj<-paste("ibd",rships[rr],sep=".")
      if(!is.null(x[[ibd.obj]])) {
        pp[[rr]] <- flagPairs(dat,x[[ibd.obj]],
                              rships.full[rr],col=rcols[rr],ellipse.coverage,
                              x[["ibd.ur"]], coverage.u, ... )
      }
    }
    # Now rbind elements of the list pp together. Use the Reduce function 
    # to succesively rbind elements.
    pp<-Reduce(rbind,pp)
  } 
  return(pp)
}

pairIDs<-function(nsub,subids=NULL) {
  # Has to pick out upper-triangular matrix column by column
  # because that's how the IBD estimates are extracted from
  # a matrix of IBD estimates in IBDest.study.
  idsub1=NULL
  idsub2=NULL
  for (i in 1:(nsub-1)){
      idsub1=c(idsub1,1:i)
      idsub2=c(idsub2,rep(i+1,i))
      #idsub1=c(idsub1,rep(i,nsub-i))
      #idsub2=c(idsub2,(i+1):nsub)
  }
  if (!is.null(subids)){
      idsub1=subids[idsub1]
      idsub2=subids[idsub2]
  }
  return(cbind(idsub1,idsub2)) 
}


pausePlotting<-function() {
  old.devask<-options("device.ask.default")
  if(old.devask==FALSE) { 
    # pause to give user time to see plot of all pairs
    devAskNewPage(ask=TRUE) 
    # Hack: draw a blank plot here -- it will be erased immediately
    plot(.5,.5,type="n")
    devAskNewPage(ask=FALSE) 
    } 
}

flagPairs<-function(dat,simdat,rship,col,coverage,simdat.unrel,coverage.u, ...) {
  # -dat is a data frame containing information on study subjects,
  #  with columns named
  #      member1: ID of the first member of the pair
  #      member2: ID of the second member of the pair
  #      pz0: estimated probability of 0 IBD for each pair
  #      pz1: estimated probability of 1 IBD for each pair
  # -simdat is the estimated iBD probs for simulated pairs
  # -rship is the relationship between simulated pairs
  # -col is the colour to plot the simulated data in 
  # -simdat.unrel is estimated IBD probs for simulated unrelated prs

  # augment the local copy of dat (local to this function) with the 
  # simulated relationship
  dat$relationship<-rship
  # create labels to appear on the plot for points identified by user
  labs<-paste(dat$member1,dat$member2,sep=":")
  # Make these labels the rownames of dat to make it easier to subset 
  #dat by the labels returned by showLabels (see below).
  rownames(dat)<-labs

  simdat<-simdat[,1:2] ; simdat.unrel<-simdat.unrel[,1:2]
  studydat<-cbind(dat$pz0,dat$pz1)

  # Infer the x- and y-limits for the plot from a prediction ellipse 
  # for the simulated pairs first
  ee<-get.ellipse(simdat,coverage=coverage)
  ee.unrel<-get.ellipse(simdat.unrel,coverage=coverage.u)
  plot(studydat[,1],studydat[,2],xlab="pz0",ylab="pz1",xlim=get.lims(ee[,1]),ylim=get.lims(ee[,2]),...)
  title(rship)
  # draw the ellipse
  lines(ee,col=col,...)
  lines(ee.unrel,col="magenta",...)
  # If simulated data is not for unrelated pairs, declare points inside 
  # the ellipse, but outside the ellipse from unrelateds, to be noteworthy
  noteworthy<-NULL
  ind<-in.ellipse(studydat,apply(simdat,2,mean),var(simdat),coverage=coverage) &
     !in.ellipse(studydat,apply(simdat.unrel,2,mean),var(simdat.unrel),coverage=coverage.u)
  noteworthy<-row.names(dat)[ind]

  #call showLabels to label interesting points and return their indices
  # NB: ** The user has to right-click the plot region in order to 
  #        exit showLabels. Warn the user of this. On Windows, may 
  #      have to flush the R-console buffer after cat() in order to get  
  #      the message to print right away. **
  cat(" ** Left-mouse-click points of interest.\n")
  cat("    Right-mouse-click the plotting region when finished **\n\n")
  if(Sys.info()["sysname"]=="Windows") flush.console()
  noteworthy<-c(noteworthy,showLabels(studydat[,1],studydat[,2],labs))
  noteworthy<-unique(noteworthy)
 
  # If you have clicked any points on the graph that are actually 
  # several points plotted on top of each other, showLabels only 
  # returns the index of the first one it finds. Call a utility  
  # function (see below) to add the rest, if they exist.
  if(length(noteworthy)>0) {
    noteworthy<-addDups(noteworthy,dat)
    out<-dat[noteworthy,]
    #Stip off labels as rownames before returning
    rownames(out)<-NULL
    return(out)
  } else {
    return(NULL)
  }

}

# Worker functions for flagPairs


get.ellipse<-function(x,coverage,...) {
  if(is.null(x)) { return(NULL) }
  if(sd(x[,1])==0) {
    # FIX ME: a bit of a kludge right now. 
    return(ellipse(0,
              scale=c(.Machine$double.eps,sd(x[,2])),
              centre=c(mean(x[,1]),mean(x[,2])),level=coverage))
  } else {
    return(ellipse(cor(x[,1],x[,2]),
              scale=c(sd(x[,1]),sd(x[,2])),
              centre=c(mean(x[,1]),mean(x[,2])),level=coverage))
  }
}

get.lims<-function(x,mult.fact=3,add.fact=0.05) {
  xr<-range(x)
  xmid<-(xr[1]+xr[2])/2
  xx<-xr[2]-xmid
  # wid<-min(mult.fact*xx,xx+add.fact)
  wid<-xx+add.fact
  return(c(xmid-wid,xmid+wid))
}

coverage.unrel<-function(coverage,npairs) {
  return(1 - (1-coverage)/npairs)
}

in.ellipse<-function(x,mu,sig,coverage) {

  # For each pair, represented by a row of the two-column matrix x, 
  # want to know if x[i,] is in the coverage*100th percentile ellipse of 
  # a normal with centre mu and variance sig; i.e., if 
  # the quadratic form Q_i = x[i,] %*% sig^{-1} %*% t(x[i,]) is less than the
  # coverage*100th percentile of a chi-squared with 2 d.f.. 
  # Can do x[i,] %*% solve(sig) for all i as tem = x %*% solve(sig). Then 
  # Q for all subjects is tem[,1]*x[,1] + tem[,2]*x[,2]

  # Have to be aware that variances could be 0. For now we replace 0's by 
  # a small number (machine precision .Machine$double.eps)
  if(sig[1,1]==0) sig[1,1]<-.Machine$double.eps
  if(sig[2,2]==0) sig[2,2]<-.Machine$double.eps
  x<-cbind(x[,1]-mu[1],x[,2]-mu[2])
  tem<-x%*%solve(sig)
  qforms<-tem[,1]*x[,1] + tem[,2]*x[,2]
  return(qforms <= qchisq(coverage,2))
}


addDups<-function(nn,dat) {
  outvec<-NULL
  allrows<-rownames(dat)
  for(i in 1:length(nn)) {
    ind<-dat[,"pz0"]==dat[nn[i],"pz0"] & dat[,"pz1"]==dat[nn[i],"pz1"]
    outvec<-c(outvec,allrows[ind])
  }
  return(outvec)
}


##function for the first plot of all study pairs clickable when 
##simulate=FALSE and non clickable when simulate=TRUE

flagAllPairs<-function(dat,simulate,ibd.ur,coverage.u,...) {
  
  labs<-paste(dat$member1,dat$member2,sep=":")
  rownames(dat)<-labs
  # plot(dat$pz0,dat$pz1,xlab="pz0", ylab="pz1",xlim=c(0,1.1),ylim=c(0,1.1),...) 
  plot(dat$pz0,dat$pz1,xlab="pz0", ylab="pz1",...) 

  if(simulate==FALSE){ # plot is clickable
    cat(" ** Left-mouse-click points of interest.\n")
    cat("    Right-mouse-click the plotting region when finished **\n\n")

   if(Sys.info()["sysname"]=="Windows") flush.console()
   noteworthy<-showLabels(dat$pz0,dat$pz1,labs)
   if(length(noteworthy)>0) {
     noteworthy<-addDups(noteworthy,dat)
     out<-dat[noteworthy,]
     #Stip off labels as rownames before returning
     rownames(out)<-NULL
     return(out)
   } else {
     return(NULL)
   }
   } else { # plot not clickable
     ee.unrel<-get.ellipse(ibd.ur[,1:2],coverage=coverage.u)
     lines(ee.unrel,col="magenta")
     return(NULL)
   }
}
