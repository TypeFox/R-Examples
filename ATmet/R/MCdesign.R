MCdesign<- function(N,k,distrib,distrib.pars,x)
{
  
  
  #MST criterion allows to compare designs
    
  #testing if distrib is missing ou incomplete
  if(missing(distrib) || length(distrib)!=k)
  {
    stop("Argument distrib is missing or not correctly specified")
  }
  
  #testing the names in distrib against the permitted names in distrib.names
  distrib.names<- list("unif","norm","t.scaled","t", "triang", "tnorm")
  if(any(distrib %in% distrib.names) ==FALSE )
  {
    stop("At least one distribution that you specify is not handled by the software, see Details...")
  }
  
  #testing the names of the variables
  if(missing(x) || length(x)!=k)
  {
    MC.names<- paste("X",1:k)
    warning("Names of variables have been replaced by 'Xn', see Details. If not relevant, please check names in x.")
  }
  else{
    MC.names<- abbreviate(x,15)
  }
  
  #######################################################
  #Checking consistency of parameters with distributions#
  #######################################################
  
  if(!missing(distrib.pars)) {
    #For gaussian distribution
    if(any(distrib=="norm") == TRUE)
    {
      ind<- which(distrib=="norm")
      for(i in ind)
      {
        if(length(distrib.pars[[i]])!=2) 
        {
          stop("incorrect number of parameters for gaussian distributions, see Details")
        }
        if(distrib.pars[[i]][2]<=0)#testing negative variance
        {
          stop("variance parameter of the Gaussian distribution should be non negative")
        }
      }
    }
  
    #For uniform distribution
    if(any(distrib=="unif") == TRUE)
    {
      ind<- which(distrib=="unif")
      for(i in ind)
      {
        if(length(distrib.pars[[i]])!=2) 
        {
          stop("incorrect number of parameters for uniform distributions, see Details")
        }
        if(as.numeric(distrib.pars[[i]][1]) > as.numeric(distrib.pars[[i]][2]))
        {
          stop("incorrect specification of boundaries for uniform distribution")
        }

      }
    }
    
    #For t distribution
    if(any(distrib=="t") == TRUE)
    {
      ind<- which(distrib=="t")
      for(i in ind)
      {
        if(length(distrib.pars[[i]])!=1) 
        {
          stop("incorrect number of parameters for t distributions, see Details")
        }
        if(distrib.pars[[i]][1]<=0)
        {
          stop("incorrect specification of parameters for t distribution")
        }
      }
    }
    
    #For t scaled distribution
    if(any(distrib=="t.scaled") == TRUE)
    {
      ind<- which(distrib=="t.scaled")
      for(i in ind)
      {
        if(length(distrib.pars[[i]])!=3) 
        {
          stop("incorrect number of parameters for t scaled distributions, see Details")
        }
        if(distrib.pars[[i]][1]<=0 || distrib.pars[[i]][3]<=0)
        {
          stop("incorrect specification of parameters for t scaled distribution")
        }
      }
    }
    
    #For triangle distribution
    if(any(distrib=="triang") == TRUE)
    {
      ind<- which(distrib=="triang")
      for(i in ind)
      {
        if(length(distrib.pars[[i]])!=3) 
        {
          stop("incorrect number of parameters for triangular distributions, see Details")
        }
      }
    }
    
    
    #For tnorm distribution
    if(any(distrib=="tnorm") == TRUE)
    {
      ind<- which(distrib=="tnorm")
      for(i in ind)
      {
        if(length(distrib.pars[[i]])!=4) 
        {
          stop("incorrect number of parameters for truncated normal distributions, see Details")
        }
      }
    }
    
    
  }
  
  #######################################################
  #End checking consistency parameters/distributions######
  #######################################################
  
  #######################################################
  #Initializing a MC design####################
  #######################################################
  
  opt.design.phys<- data.frame(matrix(0,nrow=N,ncol=k))
  
  #######################################################
  #End initializing a MC design#################
  #######################################################
  
  #######################################################
  #Creating a MC design###############################
  #######################################################
  
  if(missing(distrib.pars) || length(distrib.pars)!=k)
  {
    distrib.pars<- NULL
    warning("argument distrib.pars is missing or misspecified, default values of parameters are used instead for all the variables see Details.
            Please check distrib.pars if not relevant.")
  
    for(i in 1:k)
    {
    if(distrib[i]=="norm") opt.design.phys[,i]<- rnorm(N,0,1)
    if(distrib[i]=="unif") opt.design.phys[,i]<- runif(N,0,1)
    if(distrib[i]=="t") opt.design.phys[,i]<- rt(N,100,ncp=0)
    if(distrib[i]=="t.scaled") opt.design.phys[,i]<-rt.scaled(N, 100,0,1,ncp=0)
    if(distrib[i]=="triang") opt.design.phys[,i]<- rtri(N, min=0, max=1, mode = (min + max)/2)
    if(distrib[i]=="tnorm") opt.design.phys[,i]<- rtnorm(N, mean=0, sd=1,lower=0, upper=Inf)
  
    }
  }
  
  else{
    for(i in 1:k)
    {
      par<- as.numeric(distrib.pars[[i]])
      if(distrib[i]=="norm") opt.design.phys[,i]<- rnorm(N,par[1],par[2])
      if(distrib[i]=="unif") opt.design.phys[,i]<- runif(N,par[1],par[2])
      if(distrib[i]=="t") opt.design.phys[,i]<- rt(N,par,ncp=0)
      if(distrib[i]=="t.scaled") opt.design.phys[,i]<-rt.scaled(N, par[1],par[2],par[3],ncp=0)
      if(distrib[i]=="triang") opt.design.phys[,i]<- rtri(N, min=par[1], max=par[2], mode = par[3])
      if(distrib[i]=="tnorm") opt.design.phys[,i]<- rtnorm(N, mean=par[1], sd=par[2],lower=par[3], upper=par[4])
    }
  }
  
  #######################################################
  #End Creating a MC design###########################
  #######################################################
  
  #Naming the colums of the designs
  names(opt.design.phys)<- MC.names
  
  return(design.phys=opt.design.phys)
  
  
}
