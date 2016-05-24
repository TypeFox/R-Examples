#  File R/summary.statistics.egodata.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2016 Statnet Commons
#######################################################################
summary.statistics.egodata <- function(object,..., basis=NULL, individual=FALSE, scaleto=NULL){
  egodata <-
    if(!is.null(basis)) basis
    else get(as.character(object[[2]]), envir=environment(object))

  scaling.stats <- NULL
  scaling.pos <- c(0)
  nonscaling.stats <- c()
  nonscaling.pos <- c(0)
  
  
  for(trm in term.list.formula(object[[length(object)]])){
    if(is.call(trm)){
      init.call <- list(as.name(paste("EgoStat.", trm[[1]],sep="")),egodata=egodata)
      init.call <- c(init.call,as.list(trm[-1]))
    }else{
      init.call <- list(as.name(paste("EgoStat.", trm,sep="")),egodata=egodata)
    }
    stat<-eval(as.call(init.call), environment(object))
    if(isTRUE(attr(stat, "nonscaling"))){
      if(individual) stop("Nonscaling statistic detected. Individual contributions are meaningless.")
      nonscaling.stats <- c(nonscaling.stats, stat)
      nonscaling.pos <- c(nonscaling.pos, max(scaling.pos,nonscaling.pos) + seq_len(length(stat)))
    }else{
      scaling.stats<-cbind(scaling.stats,stat)
      scaling.pos <- c(scaling.pos, max(scaling.pos,nonscaling.pos) + seq_len(ncol(stat)))
    }
  }
  
  if(!individual){
    if(length(scaling.stats)){
      scaleto <- if(is.null(scaleto)) nrow(egodata$egos) else scaleto
      scaling.stats <- colSums(scaling.stats*egodata$egoWt)/sum(rep(egodata$egoWt,length.out=nrow(scaling.stats)))
      scaling.stats <- scaling.stats*scaleto
    }
      
    stats <- numeric(max(scaling.pos,nonscaling.pos))
    scaling.pos <- scaling.pos[scaling.pos>0]
    nonscaling.pos <- nonscaling.pos[nonscaling.pos>0]

    stats[scaling.pos] <- scaling.stats
    stats[nonscaling.pos] <- nonscaling.stats
    
    names(stats)[scaling.pos] <- names(scaling.stats)
    names(stats)[nonscaling.pos] <- names(nonscaling.stats)
    
    stats
  }else{
    rownames(scaling.stats) <- egodata$egos[[egodata$egoIDcol]]
    scaling.stats
  }
}
