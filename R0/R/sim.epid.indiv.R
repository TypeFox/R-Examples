# Name   : sim.epid.indiv
# Desc   : Individual-based model used to simulate epidemic 
# Date   : 2012/10/20 (last edit 2015/05/21)
# Author : Boelle, Obadia
###############################################################################


# Function declaration

sim.epid.indiv <- function#Influenza-like illness simulation (individual-based model)
### Generates several epidemic curves on a individual-based model
(
  beta, ##<< Contact rate in the SEIR model.
  Tmax, ##<< Maximum length of the epidemic (cases infected after this length will be truncated).
  n=1, ##<< Number of epidemics to be simulated (default is 1)
  family="poisson", ##<< Distribution of offspring (default is "poisson").
  negbin.size=NULL ##<< If family is set to "negbin", sets the size parameter of the negative binomial distribution.
) 
  
{
  
  ##note<< This is not the final version. This is the exact function as used in the manuscript (Obadia et al., 2012).
  ##It will be properly implemented to conform with other objects of the package in future releases.
  ##
  ##The epidemic is simulated using a branching process, with infinite number of susceptibles to allow for exponential growth. 
  ##The model used follows the Crump-Mode-Jagers description, with S/E/I/R description of the natural history.
  ##Latent and infectious period follow parametrized Gamma distributions typical of influenza. 
  ##An index case is first introduced, and offspring is sampled from a negative binomial distribution, with mean \eqn{beta*I} and variance \eqn{negbin.size*beta*I}, to allow for overdispersion.
  
  #epidemic matrix. Each column is an epidemic
  epid.matrix = matrix(data=0, nrow=Tmax, ncol=n)
  
  #Individual-based data.frame must be large enough to fit all cases
  if (beta <= 1.5) {
    Lmax <- 5000
  }
  else if (beta <= 2 & beta > 1.5) {
    Lmax <- 50000
  }
  else if (beta <= 3 & beta > 2) {
    Lmax <- 50000
  }
  else if (beta >3) {
    Lmax <- 75000
  }
  
  
  epid <- data.frame(matrix(data=0, nrow=Lmax, ncol=5, dimnames=list(c(1:Lmax), c("tinf","nsec","dinf","dlat","runif"))))
  
  for (i in 1:n) {
    #Memory optimization: dinf and dlat are all computed at the beginning
    epid[1:Lmax,"dinf"] <- rgamma(Lmax,shape=1/0.9**2, scale=0.9**2)
    epid[1:Lmax,"dlat"] <- rgamma(Lmax,shape=1.7**2/0.3**2, scale=0.3**2/1.7)
    
    #runif column = Time of infection of the offspring during dinf_parent
    epid[1:Lmax,"runif"] <- runif(Lmax)
    
    # add offspring
    if (family=="poisson") {
      epid[1:Lmax,"nsec"] = rpois(Lmax, lambda=epid[1:Lmax,"dinf"]*beta)
    }
    else if (family=="negbin") {
      epid[1:Lmax,"nsec"] = rnbinom(Lmax, size=negbin.size, mu=epid[1:Lmax,"dinf"]*beta)
    }
    while (epid[1,"nsec"]==0) epid[1,"nsec"] = rpois(1, lambda=epid[1,"dinf"]*beta)
    
    idx.in.epid=1
    idx.max =2
    
    while (idx.in.epid < idx.max) {
      
      #ignore infectors where tinf > Tmax
      if (epid[idx.in.epid,"tinf"] < Tmax) {
        
        #nb of offspring
        nsec = epid[idx.in.epid,"nsec"]
        if (nsec >0) {
          
          #tinf_descendant = tinf_parent + latence_parent + runif(1)*dinf_parent
          epid[idx.max:(idx.max+nsec-1),"tinf"] = epid[idx.in.epid,"tinf"] + 
            epid[idx.in.epid,"dlat"] + 
            epid[idx.max:(idx.max+nsec-1),"runif"]*epid[idx.in.epid,"dinf"]
          
          #offpsring of next infector: those after the offspring of idx.in.epid + nsec
          idx.max = idx.max + nsec
        }
      }
      #cat("idx.in.epid=",idx.in.epid,"\n idx.max=",idx.max,"\n")
      idx.in.epid = idx.in.epid+1
    }
    epid.matrix[,i] <- table(cut(epid[1:(idx.max-1),"tinf"], seq(0,Tmax,1), include.lowest=TRUE))
  }
  return(epid.matrix)
  ### A matrix with epidemics stored as columns (incidence count)
}
