corc <-
  function(dataframe,varnames,subsampsize,nbsafe=5,mixties=FALSE)
  {
    obs=dataframe[varnames]
    obs=subset(obs,apply(is.na(obs),1,sum)==0)
    nnm=length(obs[,1])
    dimension=length(varnames)
    obs=as.numeric(unlist(obs))
    u=runif(1)*(2^31-1)
    nboot=nbsafe*subsampsize^dimension
    cop=corc0( obs, nnm, dimension, subsampsize, nboot, u, mixties )
    nbootreel=cop[subsampsize^dimension+2]
    ties=cop[subsampsize^dimension+1]
    cop=cop[1:subsampsize^dimension]
    cop=array(cop, rep(subsampsize,dimension))
    cop=aperm(cop,dimension:1)
    cop=cop/((nbootreel-ties)*subsampsize)
    return(list(cop=cop,ties=ties,nsubsampreal=nbootreel,varnames=varnames,nnm=nnm))
  }
