# Utilities for dealing with MCMC output produced by *.MCMC.C functions.

ERGMM.PAR_VAR_NAMES<-c("beta","Z","sender","receiver","sociality",
                       "Z.var","Z.mean","Z.K",
                       "sender.var","receiver.var","sociality.var",
                       "dispersion")
ERGMM.PAR_LLK_NAMES<-c("beta","Z","sender","receiver","sociality","dispersion")

del.iteration<-function(mcmcsample,i){
  for(name in names(mcmcsample)){
    if(length(mcmcsample[[name]])>0){
      if(length(dim(mcmcsample[[name]]))<=1) mcmcsample[[name]]<-mcmcsample[[name]][-i]
      else if(length(dim(mcmcsample[[name]]))==2) mcmcsample[[name]]<-mcmcsample[[name]][-i,,drop=FALSE]
      else if(length(dim(mcmcsample[[name]]))==3) mcmcsample[[name]]<-mcmcsample[[name]][-i,,,drop=FALSE]
    }
  }
  mcmcsample
}

seldrop<-function(x,i){
  array(c(x),dim=dim(x)[-i])
}

as.ergmm.par.list<-function(x,...){
  class(x)<-"ergmm.par.list"
  x
}

"[.ergmm.par.list"<-function(x,i){
  if(!is.numeric(i)) stop("Index vector to the '[' operator of an ergmm.par.list must be of mode numeric or integer.")
  l<-list()
  
  for(name in names(x)){
    if(length(x[[name]])>0){
      d<-dim(x[[name]])
      if(length(d)<=1) l[[name]]<-x[[name]][i]
      else if(length(d)==2) l[[name]]<-x[[name]][i,,drop=FALSE]
      else if(length(d)==3) l[[name]]<-x[[name]][i,,,drop=FALSE]
    }
  }
  class(l)<-"ergmm.par.list"
  l
}

length.ergmm.par.list<-function(x){
  if(is.null(dim(x[[names(x)[1]]]))) length(x[[names(x)[1]]])
  else dim(x[[names(x)[1]]])[1]
}

`[[.ergmm.par.list`<-`$.ergmm.par.list`<-function(x,i){
  ## Delete its class, to keep it from recursing.
  tmp<-class(x)
  class(x)<-NULL
  if(class(i)=="character"){ ## If the index is a character, return all the draws for the corresponding variable.
    xi<-x[i][[1]]
    class(x)<-tmp
    return(xi)
  }
  else{  ## If it's a number, return a configuration with that iteration number.
    l<-list()
    ## Do NOT seldrop 1D parameters. Those actually need to become vectors.
    ## (As opposed to becoming 1*n matrices.)
    
    for(name in names(x)){
      if(length(x[[name]])>0){
        if(length(dim(x[[name]]))<=1) l[[name]]<-x[[name]][i]
        else if(length(dim(x[[name]]))==2) l[[name]]<-x[[name]][i,]
        else if(length(dim(x[[name]]))==3) l[[name]]<-seldrop(x[[name]][i,,,drop=FALSE],1)
      }
    }
    class(x)<-tmp
    return(l)
  }
}

.stack.ergmm.par.list.list<-function(x,...){
  extraneous.argcheck(...)
  mcmcsample<-list()

  for(name in names(x[[1]]))
    mcmcsample[[name]]<-abind::abind(sapply(1:length(x),
                                            function(i) x[[i]][[name]],
                                            simplify=FALSE),along=1)

  attr(mcmcsample,"breaks")<-cumsum(c(sapply(1:(length(x)
                                                 ),
                                              function(i) length(x[[i]]),
                                              simplify=FALSE)))
  class(mcmcsample)<-"ergmm.par.list"
  mcmcsample
}

unstack.ergmm.par.list<-function(x,...){
  extraneous.argcheck(...)
  mcmcList<-list()

  if(is.null(attr(x,"breaks"))){
    mcmcList[[1]]<-x
  }
  else{  
    breaks<-c(0,attr(x,"breaks"))
    
    for(i in 1:length(breaks[-1])){
      mcmcList[[i]]<-x[(breaks[i]+1):breaks[i+1]]
    }
  }
  mcmcList
}

as.mcmc.list.ergmm.par.list<-function(x,which.vars,start=1,thin=1,...){
  x<-unstack(x)
  m.l<-list()
  for(thread in 1:length(x)){
    S<-length(x[[thread]])
    m<-matrix(numeric(0),S,0)
    for(name in names(which.vars)){
      if(is.null(x[[thread]][[name]])) next
      if(length(dim(x[[thread]][[name]]))<=1){
        m2<-cbind(x[[thread]][[name]])
        colnames(m2)<-name
        m<-cbind(m,m2)
      }else if(length(dim(x[[thread]][[name]]))==2){
        m2<-x[[thread]][[name]][,c(which.vars[[name]]),drop=FALSE]
        colnames(m2)<-paste(name,sapply(which.vars[[name]],function(x) paste(x,sep='.')),sep='.')
        m<-cbind(m,m2)
      }else if(length(dim(x[[thread]][[name]]))==3){
        for(i in 1:dim(which.vars[[name]])[1]){
          i<-which.vars[[name]][i,]
          m2<-cbind(x[[thread]][[name]][,i[1],i[2]])
          colnames(m2)<-paste(name,paste(i,sep='.',collapse='.'),sep='.')
          m<-cbind(m,m2)
        }
      }
    }
    m<-mcmc(m,start=start,thin=thin)
    m.l[[thread]]<-m
  }
  eval(as.call(c(mcmc.list,m.l)))
}
