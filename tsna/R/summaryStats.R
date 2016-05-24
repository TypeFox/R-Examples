# function to apply ergm's summary formula at multiple time points
tErgmStats<-function(nd, formula, start, end, time.interval=1, aggregate.dur, rule){
  
  if(!is.networkDynamic(nd)){
    stop("the first argument to tErgmStats must be a object of class 'networkDynamic'")
  }
  
  if(is.hyper(nd)){
    stop("tErgmStats is not appropriate for hypergraphic networks")
  }
  
  if(is.multiplex(nd)){
    stop("tErgmStats is not appropriate for multiplex networks")
  }
    
    if(missing(start) | missing(end)){
      times <- get.change.times(nd)
      if (length(times) == 0) {
        warning("network does not appear to have any dynamic information. Using start=0 end=1")
        start = 0
        end = 0
      }
      times[times == Inf] <- NA
      times[times == -Inf] <- NA
      if(missing(start)){
        start <- min(times, na.rm = T)
      }
      if(missing(end)){
        end <- max(times, na.rm = T)
      }
    }
    
    # figure out the times where we will do evaluations
    times<-seq(from = start, to=end,by = time.interval)
    
    # if the formula doesn't have a '~' in front, add it
    if (!grepl('^~',formula)){
      formula<-paste('~',formula,sep='')
    }
    if(missing(aggregate.dur)){
      aggregate.dur<-0
    }
    if (missing(rule)){
      rule<-'latest'  # so that we won't have warnings about attribute processing
    }  
    
    # rquires that ergm is loaded
    # if(requireNamespace('ergm',quietly=TRUE)){
    # requireNamespace line above is prefered by CRAN, but if ergm is not loaded the necessary terms are also not loaded
    # TODO: remove ergm depends and reactivate requireNamspace code
#    if(require('ergm',quietly=TRUE)){
      stats<-lapply(times,function(t){
        # check if we are collapsing at point or interval
        if(aggregate.dur==0){
          # collapse at a point
          net<-network.collapse(nd,at=t)
        } else {
          # collapse over an interval  (slower)
          net<-network.collapse(nd,onset=t,length=aggregate.dur)
        }
        ergm::summary.statistics.formula(as.formula(paste('net',formula)))
      })
      
#    } else {
#      stop(" the ergm package could not be loaded to provide summary functions and terms")
#    }
    # rearrange list into matrix
    stats<-do.call(rbind,stats)
    return(ts(stats,start=start,end=times[length(times)],deltat=time.interval))
}

# function to provide a wrapper for calling sna measures
tSnaStats<-function(nd, snafun,start, end, time.interval=1, aggregate.dur=0, rule='latest',...){
  
  if(!is.networkDynamic(nd)){
    stop("the first argument to tSnaStats must be a object of class 'networkDynamic'")
  }
  
  # table of supported sna functions and key terms 
  # (i.e is directedness arg named 'mode' or 'gmode' )
                      # fun name  directed, type, diag
  funTerms<- matrix(c('closeness','gmode', 'VLI', 'diag',
                      'betweenness','gmode','VLI','diag',
                      'bonpow',    'gmode','VLI','diag',
                      'components', '', 'GLI', '',
                      'degree','gmode','VLI','diag',
                      'triad.census','mode','other','',
                      'connectedness','','GLI','',
                      'dyad.census','','other','',
                      'efficiency','','GLI','diag',
                      'evcent','gmode','VLI','diag',
                      'flowbet','gmode','VLI','diag',
                      'gden'     ,'mode',  'GLI', 'diag',
                      'graphcent', 'gmode','VLI','diag',
                      'grecip','','GLI','',
                      'gtrans','mode','GLI','diag',
                      'infocent','gmode','VLI','diag',
                      'hierarchy','','GLI','',
                      'loadcent','gmode','VLI','diag',
                      'lubness','','GLI','',
                      'mutuality','','GLI','',
                      'prestige','gmode','VLI','diag',
                      'centralization','mode','GLI','diag'
                     ),ncol=4,byrow=TRUE)
  
  if (!is.character(snafun)){
    stop('the "snafun" argument must be a character string giving the name of one of the supported sna package descriptive statistics')
  }
  if (!snafun%in%funTerms[,1]){
    stop('the function "', snafun,'" is not one of the sna package descriptive statistics currently supported')
  }
  
  if(missing(start) | missing(end)){
    times <- get.change.times(nd)
    if (length(times) == 0) {
      warning("network does not appear to have any dynamic information. Using start=0 end=1")
      start = 0
      end = 0
    }
    times[times == Inf] <- NA
    times[times == -Inf] <- NA
    if(missing(start)){
      start <- min(times, na.rm = T)
    }
    if(missing(end)){
      end <- max(times, na.rm = T)
    }
    
  }
  
  # figure out the times where we will do evaluations
  times<-seq(from = start, to=end,by = time.interval)
  
  args<-list(...)
  
  # rquires that sna package is loaded
  if(requireNamespace('sna',quietly=TRUE)){
    stats<-lapply(times,function(t){
      # extract the network for the time
      if(aggregate.dur==0){
        # extract at a time point (slightly faster)
        net<-network.collapse(nd,at=t)
      } else {
        # extract over an interval
        net<-network.collapse(nd,onset=t,length=aggregate.dur)
      }
      # sna functions can't handle zero-order networks
      if(network.size(net)==0){
        return(NA)
      }
      args<-c(dat=list(net),args)
      # construct appropriate args list
      # a bit messy because sometimes named mode, sometimes gmode
      directTerm<-funTerms[which(funTerms[,1]==snafun),2]
      if(!directTerm%in%names(args)){
        if(directTerm=='mode'){
          if(is.directed(net)){
            args<-c(args,mode='digraph')
          } else {
            args<-c(args,mode='graph')
          }
        } else if (directTerm=='gmode'){
          if(is.directed(net)){
            args<-c(args,gmode='digraph')
          } else {
            args<-c(args,gmode='graph')
          }
        }
        # otherwise don't add a term for directedness 
      }
      diagTerm<-funTerms[which(funTerms[,1]==snafun),4]
      if(diagTerm=='diag' && !'diag'%in%names(args)){
        if(has.loops(net)){
          args<-c(args,diag=TRUE)
        } else {
          args<-c(args,diag=FALSE)
        }
      }
      # prepend the 'sna::' namespace in case the library is not actually loated
      fun<-getExportedValue(ns = 'sna',name = snafun)
      do.call(fun,args = args)

    })
    
  } else {
    stop(" the sna package could not be loaded to provide summary functions")
  }
  # rearrange list into matrix
  # TODO: this may not work if the sizesof networks vary
  stats<-do.call(rbind,stats)
  # if it is producing one statistic per vertex, name columns appropriately
  if (funTerms[which(funTerms[,1]==snafun),3]=='VLI'){
    stats<-ts(stats,start=start,end=times[length(times)],deltat=time.interval,names=network.vertex.names(nd))
  } else {
    stats<-ts(stats,start=start,end=times[length(times)],deltat=time.interval)
  }
  return(stats)
}