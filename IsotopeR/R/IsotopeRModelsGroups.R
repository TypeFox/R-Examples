##Isotope Mixing ,Model with Measurement Error, Source Correlation, Concentration Error and residual error
##doesnt estimate third source concentration, inputs mean and variance
##outputs C:N ratios
##Jake Ferguson updated 5/4/2011
##based on original code from moore and semmens 2009


IsotopeRfullgroup <- "model {

  ##################################
  ##estimate the measurement error##
  ##################################
  for(iso1 in 2:num.iso) {
    for(iso2 in 1:(iso1-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso2 in 2:num.iso) {
    for(iso1 in 1:(iso2-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso in 1:num.iso) {
    tauZ[iso,iso] ~ dgamma(.001,.001)
  }

  cov.ME <- inverse(tauZ);

  for(i in 1:Nz) {
    Z[i,] ~ dmnorm(muz,tauZ);
  }
    
  ###############################
  ##estimate the concentrations##
  ###############################
  for(source in 1:(num.sources)) {
  ##covariance matrix
  for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        D.tau[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        D.tau[sourcex,sourcey,source] <- 0
      }
    }    
    for(source2 in 1:num.iso) {
      D.tau[source2,source2,source] ~ dexp(0.001) #dgamma(0.001, 0.001)
    }

	##draws subsource means and fits to data
#     for(sub in 1:subcd.vec[source]) { 
# 		subcd[source, sub, 1:num.iso] ~ dmnorm(dmu.prior.mu[,source], dmu.prior.tau)
# 	}
	##draws subsource means and fits to data
    for(sub in 1:subcd.vec[source]) { 
	  for(iso in 1:num.iso) {
		  subcd[source, sub, iso] ~ dunif(0,100)
		}
	}
	
	for(iso in 1:num.iso) {
		mu.conc2[iso, source] <- mean( abs(subcd[source, 1:subcd.vec[source], iso]))
	}
		
	for(sub in 1:subcd.vec[source]) {
		for(index  in subcd.samples[source,sub,1]:subcd.samples[source,sub,2]) {
			cd.mat[index, ] ~ dmnorm(subcd[source,sub, ], D.tau[ , ,source]);
		}  
  }

  }
  
  
  ###########################
  ##Proportion estimation##
  ###########################
  ## this is the global mean
  for(i in 1:num.sources) {mu[i] ~ dnorm(alpha.clr[i], 0.001)} 
  pop.invSig2 ~ dgamma(.1,.1)
  pop.var <- 1/pop.invSig2
  for(source in 1:num.sources) {
    p.transform[source] ~ dnorm(mu[source],pop.invSig2);
    p.exp.transform[source] <- exp(p.transform[source]);
  }  
  p.mean.tot <- sum(p.exp.transform[]); 
   
  ####################################
  ##this is the subpopulations means##
  ####################################
  subpop.invSig2 ~ dgamma(0.1,0.1)
  subpop.var <- 1/subpop.invSig2
  for(group in 1:num.groups) {
    for(source in 1:num.sources) {
      p.group.clr[group,source] ~ dnorm(p.transform[source], subpop.invSig2)
      p.group.transform[group,source] <- exp(p.group.clr[group,source])
    }
    p.group.tot[group] <- sum(p.group.transform[group,]); 
  }
  
  ##generate individuals draws from the group mean
  ind.invSig2 ~ dgamma(.1,.1)
  ind.var <- 1/ind.invSig2
  for(group in 1:num.groups) {
    for(i in groupnum.mat[group,1]:groupnum.mat[group,2]) { 
      for(source in 1:num.sources) {
		    p.ind[i,source] ~ dnorm(p.group.clr[group,source], ind.invSig2);
		    exp.p[i,source] <- exp(p.ind[i,source]);
      }
    }
  }
    
  ##CLR math: This does the back-transform to get back to proportions
  for(source in 1:(num.sources-1)) {
    p.pop[source] <- p.exp.transform[source]/p.mean.tot;
  }
  p.pop[num.sources] <- 1-sum(p.pop[1:(num.sources-1)]); 
  
  for(group in 1:num.groups) {
    for(source in 1:(num.sources-1)) {
      p.group[group,source] <- p.group.transform[group,source]/p.group.tot[group]
    }
    p.group[group,num.sources] <- 1-sum(p.group[group,1:(num.sources-1)]); 
  }
    
  ##individual p's
  for(i in 1:N) {
    p.tot[i] <- sum(exp.p[i,1:num.sources]);
    for(source in 1:(num.sources-1)) {
	     p[i,source] <- exp.p[i,source]/p.tot[i];
    }
    p[i,num.sources] <- 1-sum(p[i,1:(num.sources-1)]);
   
    ##rescale p.pop for concentration dependence
    p.inddenom[i,1:num.iso] <- mu.conc2%*%p[i,]
    for(iso in 1:num.iso) {
      for(source in 1:num.sources) {
		    pIso.ind[i,iso,source]	<- mu.conc2[iso,source]*p[i,source]/p.inddenom[i,iso]
      }
    }
 }#end for i


  #####################
  ##Source Estimation##
  #####################
  ##estimate sources
  for(source in 1:num.sources) { 
   ##covariance matrix
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
		tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
		tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(source2 in 1:num.iso) {
      tau.source.temp[source2,source2,source] ~ dgamma(.0010,.0010)
    }
       
    ##build source correlation matrix
   for(sourcex in 1:num.iso) {
      for(sourcey in 1:num.iso) {
      rho.source[sourcex,sourcey,source] ~ dunif(-0.99,0.99);
    }
    }
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        rho.mat[source, sourcex,sourcey] <- rho.source[sourcex, sourcey, source]*rho.flag
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        rho.mat[source,sourcex,sourcey] <- rho.mat[source, sourcey, sourcex];
      }
    }
    for(source2 in 1:num.iso) {
      rho.mat[source,source2,source2] <- 1;
    }      
    
    cov.source[1:num.iso,1:num.iso,source] <- inverse(tau.source.temp[,,source])    
    tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[source,,]%*%cov.source[,,source]) +  cov.ME  )
 
    ##draws subsource means and fit to data
    for(sub in 1:subsource.vec[source]) { 
		subsource[source, sub, 1:num.iso] ~ dmnorm(mu.prior.mu, mu.prior.cov)
	}
	
	for(iso in 1:num.iso) {
		mu.source[source, iso] <- mean( subsource[source, 1:subsource.vec[source], iso])
	}
		
	for(sub in 1:subsource.vec[source]) {
		for(index  in subsource.samples[source,sub,1]:subsource.samples[source,sub,2]) {
			source.mat[index, ] ~ dmnorm(subsource[source,sub, ], tau.source[ , ,source]);
		}  
     }

  }#end for sources
  
  
  ####################### 
  ##draw residual error##
  #######################
  for(i in 1:(num.iso-1)) {
	  for(j in (i+1):num.iso) {
		res.tau[i,j] <- 0;
		res.tau[j,i] <- 0;
		}
	}
  for(iso in 1:num.iso) {
    res.tau[iso,iso] ~ dgamma(1e-3,1e-3)#dexp(1/1000)#dunif(0,20)#dgamma(10,10)#dunif(0,20);#dexp(1);
  }
  res.err[1:num.iso,1:num.iso] <- inverse(res.tau)

  ##rescale sources by p
  for(i in 1:N) {
    ##rescale covariance and include fractionation
    for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
        covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source])*pIso.ind[i,iso,source]
      }
      for(isox in 2:(num.iso)) {
		for(isoy in 1:(isox-1)) {
		  covfrac.source[isox,isoy,source,i] <- 0
		}
      }
      for(isoy in 2:(num.iso)) {
		for(isox in 1:(isoy-1)) {
		  covfrac.source[isox,isoy,source,i] <- 0
		}
      }
  
      obscov.mat[1:num.iso,1:num.iso,source,i] <- (covfrac.source[,,source,i]%*%rho.mat[source,,]%*%covfrac.source[,,source,i] +  cov.ME + res.err)  
  
    }#end for source

    for(x in 1:num.iso) {
      for(y in 1:num.iso) {
		sumobscov.mat[x,y,i] <- sum(obscov.mat[x,y,1:num.sources,i])  
      }
    }
    for(iso in 1:num.iso) {
      mu.mix[i,iso] <- pIso.ind[i,iso,]%*%(mu.source[,iso])
    }
  }#end for i

  ##get the sd's for jack
  for(iso in 1:num.iso) {
  
    sd.res[iso] <- sqrt(res.err[iso,iso])
    sd.me[iso] <- sqrt(cov.ME[iso,iso])
            
    for(source in 1:num.sources) {
        sd.source[source,iso] <- sqrt(cov.source[iso,iso,source])
        sd.conc[source,iso] <- 1/sqrt(D.tau[iso,iso,source])
    }
  }

  mu.conc	<- t(mu.conc2)

  ##calculate the likelihoods for the N individuals.
  for(ind in 1:N) {
    mix.prcsn[1:num.iso,1:num.iso,ind] <- inverse(sumobscov.mat[,,ind]  )
    for(j in 1:ind.counts[ind]) {      
      ind.array[1:num.iso,ind,j] ~ dmnorm(mu.mix[ind,1:num.iso], mix.prcsn[1:num.iso,1:num.iso, ind]);    
    }
  }

}"


IsotopeRnoconcgroup <- "model {

  ##################################
  ##estimate the measurement error##
  ##################################
  for(iso1 in 2:num.iso) {
    for(iso2 in 1:(iso1-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso2 in 2:num.iso) {
    for(iso1 in 1:(iso2-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso in 1:num.iso) {
    tauZ[iso,iso] ~ dgamma(.001,.001)
  }

  cov.ME <- inverse(tauZ);

  for(i in 1:Nz) {
    Z[i,] ~ dmnorm(muz,tauZ);
  }
    
  ###############################
  ##estimate the concentrations##
  ###############################
  for(source in 1:(num.sources)) {
      for(iso in 1:num.iso) {
        mu.conc2[iso,source] <- 1
      }
  }
  
  
  ###########################
  ##Proportion estimamation##
  ###########################
  ## this is the global mean
  for(i in 1:num.sources) {mu[i] ~ dnorm(alpha.clr[i], 0.01);} 
  pop.invSig2 ~ dgamma(.01, .01);
  pop.var <- 1/pop.invSig2;
  for(source in 1:num.sources) {
    p.transform[source] ~ dnorm(mu[source], pop.invSig2);
  }  

	subpop.invSig2 ~ dgamma(0.1,0.1)
    for(group in 1:num.groups) {
		for(source in 1:num.sources) {
			p.group.clr[group,source] ~ dnorm(p.transform[source], subpop.invSig2)
			p.group.transform[group,source] <- exp(p.group.clr[group,source])
		}
    p.group.tot[group] <- sum(p.group.transform[group,]); 
  }

for(group in 1:num.groups) {
    for(source in 1:(num.sources-1)) {
      p.group[group,source] <- p.group.transform[group,source]/p.group.tot[group]
    }
    p.group[group,num.sources] <- 1-sum(p.group[group,1:(num.sources-1)]); 
  }

  
  ##generate individuals draws from the global mean
  ind.invSig2 ~ dgamma(.01, .01);
  ind.var <- 1/ind.invSig2;
  for(group in 1:num.groups) {
    for(i in groupnum.mat[group,1]:groupnum.mat[group,2]) { 
      for(source in 1:num.sources) {
        p.ind[i,source] ~ dnorm(p.group.clr[group,source], ind.invSig2);
        exp.p[i,source] <- exp(p.ind[i,source]);
      }
    }
  }

    
  ##CLR math: This does the back-transform to get back to proportions
  for(source in 1:num.sources) {
    p.exp.transform[source] <- exp(p.transform[source]);
  }
  p.mean.tot <- sum(p.exp.transform[]);
  
  for(source in 1:(num.sources-1)) {
    p.pop[source] <- p.exp.transform[source]/p.mean.tot;
  }
  p.pop[num.sources] <- 1-sum(p.pop[1:(num.sources-1)]); 
  
  ##rescale p.pop for concentration dependence
  for(iso in 1:num.iso) {
    for(source in 1:num.sources) {
      pIso.pop[iso,source]  <- mu.conc2[iso,source]*p.pop[source]/sum(mu.conc2[iso,]*p.pop) #p.popdenom[iso]
    }
  }
  
  ##individual p's
  for(i in 1:N) {
    p.tot[i] <- sum(exp.p[i,1:num.sources]);
      for(source in 1:(num.sources-1)) {
    p[i,source] <- exp.p[i,source]/p.tot[i];
      }
    p[i,num.sources] <- 1-sum(p[i,1:(num.sources-1)]);
   
    ##rescale p.pop for concentration dependence
    p.inddenom[i,1:num.iso] <- mu.conc2%*%p[i,]
    for(iso in 1:num.iso) {
      for(source in 1:num.sources) {
        pIso.ind[i,iso,source]  <- mu.conc2[iso,source]*p[i,source]/p.inddenom[i,iso]
      }
    }
 }#end for i
  
  #####################
  ##Source Estimation##
  #####################
  ##estimate sources
  for(source in 1:num.sources) { 
   ##covariance matrix
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(source2 in 1:num.iso) {
      tau.source.temp[source2,source2,source] ~ dgamma(.0010,.0010)
    }
       
    ##build source correlation matrix
   for(sourcex in 1:num.iso) {
      for(sourcey in 1:num.iso) {
 		  rho.source[sourcex,sourcey,source] ~ dunif(-0.99,0.99);
		}
    }
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        rho.mat[source, sourcex,sourcey] <- rho.source[sourcex, sourcey, source]*rho.flag
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        rho.mat[source,sourcex,sourcey] <- rho.mat[source, sourcey, sourcex];
      }
    }
    for(source2 in 1:num.iso) {
      rho.mat[source,source2,source2] <- 1;
    }      

    cov.source[1:num.iso,1:num.iso,source] <- inverse(tau.source.temp[,,source])
    tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[source,,]%*%cov.source[,,source]) +  cov.ME )
 
    ##draws subsource means and fit to data
    for(sub in 1:subsource.vec[source]) { 
		subsource[source, sub, 1:num.iso] ~ dmnorm(mu.prior.mu, mu.prior.cov)
	}
	
	for(iso in 1:num.iso) {
		mu.source[source, iso] <- mean( subsource[source, 1:subsource.vec[source], iso])
	}
		
	for(sub in 1:subsource.vec[source]) {
		for(index  in subsource.samples[source,sub,1]:subsource.samples[source,sub,2]) {
			source.mat[index, ] ~ dmnorm(subsource[source,sub, ], tau.source[ , ,source]);
		}  
     }

	}#end for sources
  
  
  ####################### 
  ##draw residual error##
  #######################
  for(i in 1:(num.iso-1)) {
	  for(j in (i+1):num.iso) {
		res.tau[i,j] <- 0;
		res.tau[j,i] <- 0;
		}
	}
  for(iso in 1:num.iso) {
    res.tau[iso,iso] ~ dgamma(1e-3,1e-3)#dexp(1/1000)#dunif(0,20)#dgamma(10,10)#dunif(0,20);#dexp(1);
  }
  res.err[1:num.iso,1:num.iso] <- inverse(res.tau)

  ##rescale sources by p
  for(i in 1:N) {
    ##rescale covariance and include fractionation
    for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
        covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source])*pIso.ind[i,iso,source]
      }
      for(isox in 2:(num.iso)) {
    for(isoy in 1:(isox-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
      for(isoy in 2:(num.iso)) {
    for(isox in 1:(isoy-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
  
      obscov.mat[1:num.iso,1:num.iso,source,i] <- (covfrac.source[,,source,i]%*%rho.mat[source,,]%*%covfrac.source[,,source,i] +  cov.ME + res.err)  
  
    }#end for source

    for(x in 1:num.iso) {
      for(y in 1:num.iso) {
    sumobscov.mat[x,y,i] <- sum(obscov.mat[x,y,1:num.sources,i])  
      }
    }
    for(iso in 1:num.iso) {
      mu.mix[i,iso] <- pIso.ind[i,iso,]%*%(mu.source[,iso])
    }
  }#end for i

  ##get the sd's for jack
  for(iso in 1:num.iso) {
  
    sd.res[iso] <- sqrt(res.err[iso,iso])
    sd.me[iso] <- sqrt(cov.ME[iso,iso])
            
    for(source in 1:num.sources) {

        sd.source[source,iso] <- sqrt(cov.source[iso,iso,source])
#         sd.conc[source,iso] <- 1/sqrt(D.tau[iso,iso,source])

    }
  }

  mu.conc	<- t(mu.conc2)

  ##calculate the likelihoods for the N individuals.
  for(ind in 1:N) {
    mix.prcsn[1:num.iso,1:num.iso,ind] <- inverse(sumobscov.mat[,,ind]  )
    for(j in 1:ind.counts[ind]) {      
      ind.array[1:num.iso,ind,j] ~ dmnorm(mu.mix[ind,1:num.iso], mix.prcsn[1:num.iso,1:num.iso, ind]);    
    }
  }
}"





IsotopeRnoconcgroup <- "model {

  ##################################
  ##estimate the measurement error##
  ##################################
  for(iso1 in 2:num.iso) {
    for(iso2 in 1:(iso1-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso2 in 2:num.iso) {
    for(iso1 in 1:(iso2-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso in 1:num.iso) {
    tauZ[iso,iso] ~ dgamma(.001,.001)
  }

  cov.ME <- inverse(tauZ);

  for(i in 1:Nz) {
    Z[i,] ~ dmnorm(muz,tauZ);
  }
    
  ###############################
  ##estimate the concentrations##
  ###############################
  for(source in 1:(num.sources)) {
      for(iso in 1:num.iso) {
        mu.conc2[iso,source] <- 1
      }
  }
  
  
  ###########################
  ##Proportion estimamation##
  ###########################
  ## this is the global mean
  for(i in 1:num.sources) {mu[i] ~ dnorm(alpha.clr[i], 0.01);} 
  pop.invSig2 ~ dgamma(.01, .01);
  pop.var <- 1/pop.invSig2;
  for(source in 1:num.sources) {
    p.transform[source] ~ dnorm(mu[source], pop.invSig2);
  }  

	subpop.invSig2 ~ dgamma(0.1,0.1)
    for(group in 1:num.groups) {
		for(source in 1:num.sources) {
			p.group.clr[group,source] ~ dnorm(p.transform[source], subpop.invSig2)
			p.group.transform[group,source] <- exp(p.group.clr[group,source])
		}
    p.group.tot[group] <- sum(p.group.transform[group,]); 
  }
  
  ##generate individuals draws from the global mean
  ind.invSig2 ~ dgamma(.01, .01);
  ind.var <- 1/ind.invSig2;
  for(group in 1:num.groups) {
    for(i in groupnum.mat[group,1]:groupnum.mat[group,2]) { 
      for(source in 1:num.sources) {
        p.ind[i,source] ~ dnorm(p.group.clr[group,source], ind.invSig2);
        exp.p[i,source] <- exp(p.ind[i,source]);
      }
    }
  }
    
  ##CLR math: This does the back-transform to get back to proportions
  for(source in 1:num.sources) {
    p.exp.transform[source] <- exp(p.transform[source]);
  }
  p.mean.tot <- sum(p.exp.transform[]);
  
  for(source in 1:(num.sources-1)) {
    p.pop[source] <- p.exp.transform[source]/p.mean.tot;
  }
  p.pop[num.sources] <- 1-sum(p.pop[1:(num.sources-1)]); 
  
  for(group in 1:num.groups) {
    for(source in 1:(num.sources-1)) {
      p.group[group,source] <- p.group.transform[group,source]/p.group.tot[group]
    }
    p.group[group,num.sources] <- 1-sum(p.group[group,1:(num.sources-1)]); 
  }

  ##rescale p.pop for concentration dependence
  for(iso in 1:num.iso) {
    for(source in 1:num.sources) {
      pIso.pop[iso,source]  <- mu.conc2[iso,source]*p.pop[source]/sum(mu.conc2[iso,]*p.pop) #p.popdenom[iso]
    }
  }
  
  ##individual p's
  for(i in 1:N) {
    p.tot[i] <- sum(exp.p[i,1:num.sources]);
      for(source in 1:(num.sources-1)) {
    p[i,source] <- exp.p[i,source]/p.tot[i];
      }
    p[i,num.sources] <- 1-sum(p[i,1:(num.sources-1)]);
   
    ##rescale p.pop for concentration dependence
    p.inddenom[i,1:num.iso] <- mu.conc2%*%p[i,]
    for(iso in 1:num.iso) {
      for(source in 1:num.sources) {
        pIso.ind[i,iso,source]  <- mu.conc2[iso,source]*p[i,source]/p.inddenom[i,iso]
      }
    }
 }#end for i
  
  #####################
  ##Source Estimation##
  #####################
  ##estimate sources
  for(source in 1:num.sources) { 
   ##covariance matrix
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(source2 in 1:num.iso) {
      tau.source.temp[source2,source2,source] ~ dgamma(.0010,.0010)
    }
       
    ##build source correlation matrix
    ##build source correlation matrix
   for(sourcex in 1:num.iso) {
      for(sourcey in 1:num.iso) {
 		  rho.source[sourcex,sourcey,source] ~ dunif(-0.99,0.99);
		}
    }
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        rho.mat[source, sourcex,sourcey] <- rho.source[sourcex, sourcey, source]*rho.flag
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        rho.mat[source,sourcex,sourcey] <- rho.mat[source, sourcey, sourcex];
      }
    }
    for(source2 in 1:num.iso) {
      rho.mat[source,source2,source2] <- 1;
    }      

    cov.source[1:num.iso,1:num.iso,source] <- inverse(tau.source.temp[,,source])
    tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[source,,]%*%cov.source[,,source]) +  cov.ME )
 
    ##draws subsource means and fit to data
    for(sub in 1:subsource.vec[source]) { 
		subsource[source, sub, 1:num.iso] ~ dmnorm(mu.prior.mu, mu.prior.cov)
	}
	
	for(iso in 1:num.iso) {
		mu.source[source, iso] <- mean( subsource[source, 1:subsource.vec[source], iso])
	}
		
	for(sub in 1:subsource.vec[source]) {
		for(index  in subsource.samples[source,sub,1]:subsource.samples[source,sub,2]) {
			source.mat[index, ] ~ dmnorm(subsource[source,sub, ], tau.source[ , ,source]);
		}  
     }
  }#end for sources
  
  
  ####################### 
  ##draw residual error##
  #######################
  for(i in 1:(num.iso-1)) {
	  for(j in (i+1):num.iso) {
		res.tau[i,j] <- 0;
		res.tau[j,i] <- 0;
		}
	}
  for(iso in 1:num.iso) {
    res.tau[iso,iso] ~ dgamma(1e-3,1e-3)#dexp(1/1000)#dunif(0,20)#dgamma(10,10)#dunif(0,20);#dexp(1);
  }
  res.err[1:num.iso,1:num.iso] <- inverse(res.tau)

  ##rescale sources by p
  for(i in 1:N) {
    ##rescale covariance and include fractionation
    for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
        covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source])*pIso.ind[i,iso,source]
      }
      for(isox in 2:(num.iso)) {
    for(isoy in 1:(isox-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
      for(isoy in 2:(num.iso)) {
    for(isox in 1:(isoy-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
  
      obscov.mat[1:num.iso,1:num.iso,source,i] <- (covfrac.source[,,source,i]%*%rho.mat[source,,]%*%covfrac.source[,,source,i] +  cov.ME + res.err)  
  
    }#end for source

    for(x in 1:num.iso) {
      for(y in 1:num.iso) {
    sumobscov.mat[x,y,i] <- sum(obscov.mat[x,y,1:num.sources,i])  
      }
    }
    for(iso in 1:num.iso) {
      mu.mix[i,iso] <- pIso.ind[i,iso,]%*%(mu.source[,iso])
    }
  }#end for i

  ##get the sd's for jack
  for(iso in 1:num.iso) {
  
    sd.res[iso] <- sqrt(res.err[iso,iso])
    sd.me[iso] <- sqrt(cov.ME[iso,iso])
            
    for(source in 1:num.sources) {

        sd.source[source,iso] <- sqrt(cov.source[iso,iso,source])
#         sd.conc[source,iso] <- 1/sqrt(D.tau[iso,iso,source])

    }
  }

  mu.conc	<- t(mu.conc2)

  ##calculate the likelihoods for the N individuals.
  for(ind in 1:N) {
    mix.prcsn[1:num.iso,1:num.iso,ind] <- inverse(sumobscov.mat[,,ind]  )
    for(j in 1:ind.counts[ind]) {      
      ind.array[1:num.iso,ind,j] ~ dmnorm(mu.mix[ind,1:num.iso], mix.prcsn[1:num.iso,1:num.iso, ind]);    
    }
  }

}"



IsotopeRnoconcnomegroup <- "model {
    
  ###############################
  ##estimate the concentrations##
  ###############################
  for(source in 1:(num.sources)) {
      for(iso in 1:num.iso) {
        mu.conc2[iso,source] <- 1
      }
  }
  
  
  ###########################
  ##Proportion estimamation##
  ###########################
  ## this is the global mean
  for(i in 1:num.sources) {mu[i] ~ dnorm(alpha.clr[i], 0.01);} 
  pop.invSig2 ~ dgamma(.01, .01);
  pop.var <- 1/pop.invSig2;
  for(source in 1:num.sources) {
    p.transform[source] ~ dnorm(mu[source], pop.invSig2);
  }  

	subpop.invSig2 ~ dgamma(0.1,0.1)
    for(group in 1:num.groups) {
		for(source in 1:num.sources) {
			p.group.clr[group,source] ~ dnorm(p.transform[source], subpop.invSig2)
			p.group.transform[group,source] <- exp(p.group.clr[group,source])
		}
    p.group.tot[group] <- sum(p.group.transform[group,]); 
  }
  
  ##generate individuals draws from the global mean
  ind.invSig2 ~ dgamma(.01, .01);
  ind.var <- 1/ind.invSig2;
  for(group in 1:num.groups) {
    for(i in groupnum.mat[group,1]:groupnum.mat[group,2]) { 
      for(source in 1:num.sources) {
        p.ind[i,source] ~ dnorm(p.group.clr[group,source], ind.invSig2);
        exp.p[i,source] <- exp(p.ind[i,source]);
      }
    }
  }
    
  ##CLR math: This does the back-transform to get back to proportions
  for(source in 1:num.sources) {
    p.exp.transform[source] <- exp(p.transform[source]);
  }
  p.mean.tot <- sum(p.exp.transform[]);
  
  for(source in 1:(num.sources-1)) {
    p.pop[source] <- p.exp.transform[source]/p.mean.tot;
  }
  p.pop[num.sources] <- 1-sum(p.pop[1:(num.sources-1)]); 
  
  for(group in 1:num.groups) {
    for(source in 1:(num.sources-1)) {
      p.group[group,source] <- p.group.transform[group,source]/p.group.tot[group]
    }
    p.group[group,num.sources] <- 1-sum(p.group[group,1:(num.sources-1)]); 
  }

  ##rescale p.pop for concentration dependence
  for(iso in 1:num.iso) {
    for(source in 1:num.sources) {
      pIso.pop[iso,source]  <- mu.conc2[iso,source]*p.pop[source]/sum(mu.conc2[iso,]*p.pop) #p.popdenom[iso]
    }
  }
  
  ##individual p's
  for(i in 1:N) {
    p.tot[i] <- sum(exp.p[i,1:num.sources]);
      for(source in 1:(num.sources-1)) {
    p[i,source] <- exp.p[i,source]/p.tot[i];
      }
    p[i,num.sources] <- 1-sum(p[i,1:(num.sources-1)]);
   
    ##rescale p.pop for concentration dependence
    p.inddenom[i,1:num.iso] <- mu.conc2%*%p[i,]
    for(iso in 1:num.iso) {
      for(source in 1:num.sources) {
        pIso.ind[i,iso,source]  <- mu.conc2[iso,source]*p[i,source]/p.inddenom[i,iso]
      }
    }
 }#end for i
  
  #####################
  ##Source Estimation##
  #####################
  ##estimate sources
  for(source in 1:num.sources) { 
   ##covariance matrix
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(source2 in 1:num.iso) {
      tau.source.temp[source2,source2,source] ~ dgamma(.0010,.0010)
    }
       
    ##build source correlation matrix
    ##build source correlation matrix
   for(sourcex in 1:num.iso) {
      for(sourcey in 1:num.iso) {
 		  rho.source[sourcex,sourcey,source] ~ dunif(-0.99,0.99);
		}
    }
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        rho.mat[source, sourcex,sourcey] <- rho.source[sourcex, sourcey, source]*rho.flag
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        rho.mat[source,sourcex,sourcey] <- rho.mat[source, sourcey, sourcex];
      }
    }
    for(source2 in 1:num.iso) {
      rho.mat[source,source2,source2] <- 1;
    }      

    cov.source[1:num.iso,1:num.iso,source] <- inverse(tau.source.temp[,,source])
    tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[source,,]%*%cov.source[,,source])  )
 
    ##draws subsource means and fit to data
    for(sub in 1:subsource.vec[source]) { 
		subsource[source, sub, 1:num.iso] ~ dmnorm(mu.prior.mu, mu.prior.cov)
	}
	
	for(iso in 1:num.iso) {
		mu.source[source, iso] <- mean( subsource[source, 1:subsource.vec[source], iso])
	}
		
	for(sub in 1:subsource.vec[source]) {
		for(index  in subsource.samples[source,sub,1]:subsource.samples[source,sub,2]) {
			source.mat[index, ] ~ dmnorm(subsource[source,sub, ], tau.source[ , ,source]);
		}  
     }
  }#end for sources
  
  
  ####################### 
  ##draw residual error##
  #######################
  for(i in 1:(num.iso-1)) {
	  for(j in (i+1):num.iso) {
		res.tau[i,j] <- 0;
		res.tau[j,i] <- 0;
		}
	}
  for(iso in 1:num.iso) {
    res.tau[iso,iso] ~ dgamma(1e-3,1e-3)#dexp(1/1000)#dunif(0,20)#dgamma(10,10)#dunif(0,20);#dexp(1);
  }
  res.err[1:num.iso,1:num.iso] <- inverse(res.tau)

  ##rescale sources by p
  for(i in 1:N) {
    ##rescale covariance and include fractionation
    for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
        covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source] )*pIso.ind[i,iso,source]
      }
      for(isox in 2:(num.iso)) {
    for(isoy in 1:(isox-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
      for(isoy in 2:(num.iso)) {
    for(isox in 1:(isoy-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
  
      obscov.mat[1:num.iso,1:num.iso,source,i] <- (covfrac.source[,,source,i]%*%rho.mat[source,,]%*%covfrac.source[,,source,i] +  res.err)  
  
    }#end for source

    for(x in 1:num.iso) {
      for(y in 1:num.iso) {
    sumobscov.mat[x,y,i] <- sum(obscov.mat[x,y,1:num.sources,i])  
      }
    }
    for(iso in 1:num.iso) {
      mu.mix[i,iso] <- pIso.ind[i,iso,]%*%(mu.source[,iso])
    }
  }#end for i

  ##get the sd's for jack
  for(iso in 1:num.iso) {
  
    sd.res[iso] <- sqrt(res.err[iso,iso])
#     sd.me[iso] <- sqrt(cov.ME[iso,iso])
            
    for(source in 1:num.sources) {

        sd.source[source,iso] <- sqrt(cov.source[iso,iso,source])
#         sd.conc[source,iso] <- 1/sqrt(D.tau[iso,iso,source])

    }
  }

  mu.conc	<- t(mu.conc2)

  ##calculate the likelihoods for the N individuals.
  for(ind in 1:N) {
    mix.prcsn[1:num.iso,1:num.iso,ind] <- inverse(sumobscov.mat[,,ind]  )
    for(j in 1:ind.counts[ind]) {      
      ind.array[1:num.iso,ind,j] ~ dmnorm(mu.mix[ind,1:num.iso], mix.prcsn[1:num.iso,1:num.iso, ind]);    
    }
  }

}"






IsotopeRnomegroup <- "model {
    
  ###############################
  ##estimate the concentrations##
  ###############################
  for(source in 1:(num.sources)) {
  ##covariance matrix
  for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        D.tau[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        D.tau[sourcex,sourcey,source] <- 0
      }
    }    
    for(source2 in 1:num.iso) {
      D.tau[source2,source2,source] ~ dexp(0.001)#dgamma(0.001, 0.001)
    }

	##draws subsource means and fits to data
#     for(sub in 1:subcd.vec[source]) { 
# 		subcd[source, sub, 1:num.iso] ~ dmnorm(dmu.prior.mu[,source], dmu.prior.tau)
# 	}
	##draws subsource means and fits to data
    for(sub in 1:subcd.vec[source]) { 
	  for(iso in 1:num.iso) {
		  subcd[source, sub, iso] ~ dunif(0,100)
		}
	}
	
	for(iso in 1:num.iso) {
		mu.conc2[iso, source] <- mean( abs(subcd[source, 1:subcd.vec[source], iso]))
	}
		
	for(sub in 1:subcd.vec[source]) {
		for(index  in subcd.samples[source,sub,1]:subcd.samples[source,sub,2]) {
			cd.mat[index, ] ~ dmnorm(subcd[source,sub, ], D.tau[ , ,source]);
		}  
     }

  }
  
  
  ###########################
  ##Proportion estimamation##
  ###########################
  ## this is the global mean
  for(i in 1:num.sources) {mu[i] ~ dnorm(alpha.clr[i], 0.01);} 
  pop.invSig2 ~ dgamma(.01, .01);
  pop.var <- 1/pop.invSig2;
  for(source in 1:num.sources) {
    p.transform[source] ~ dnorm(mu[source], pop.invSig2);
  }  

	subpop.invSig2 ~ dgamma(0.1,0.1)
    for(group in 1:num.groups) {
	  for(source in 1:num.sources) {
		p.group.clr[group,source] ~ dnorm(p.transform[source], subpop.invSig2)
		p.group.transform[group,source] <- exp(p.group.clr[group,source])
	  }
	  p.group.tot[group] <- sum(p.group.transform[group,]); 
  }
  
  ##generate individuals draws from the global mean
  ind.invSig2 ~ dgamma(.01, .01);
  ind.var <- 1/ind.invSig2;
  for(group in 1:num.groups) {
    for(i in groupnum.mat[group,1]:groupnum.mat[group,2]) { 
      for(source in 1:num.sources) {
		    p.ind[i,source] ~ dnorm(p.group.clr[group,source], ind.invSig2);
		    exp.p[i,source] <- exp(p.ind[i,source]);
      }
    }
  }
  
    
  ##CLR math: This does the back-transform to get back to proportions
  for(source in 1:num.sources) {
    p.exp.transform[source] <- exp(p.transform[source]);
  }
  p.mean.tot <- sum(p.exp.transform[]);
  
  for(source in 1:(num.sources-1)) {
    p.pop[source] <- p.exp.transform[source]/p.mean.tot;
  }
  p.pop[num.sources] <- 1-sum(p.pop[1:(num.sources-1)]); 
  
  for(group in 1:num.groups) {
    for(source in 1:(num.sources-1)) {
      p.group[group,source] <- p.group.transform[group,source]/p.group.tot[group]
    }
    p.group[group, num.sources] <- 1-sum(p.group[group,1:(num.sources-1)]); 
  }

  ##rescale p.pop for concentration dependence
  for(iso in 1:num.iso) {
    for(source in 1:num.sources) {
      pIso.pop[iso,source]  <- mu.conc2[iso,source]*p.pop[source]/sum(mu.conc2[iso,]*p.pop) #p.popdenom[iso]
    }
  }

  ##individual p's
  for(i in 1:N) {
    p.tot[i] <- sum(exp.p[i,1:num.sources]);
      for(source in 1:(num.sources-1)) {
        p[i,source] <- exp.p[i,source]/p.tot[i];
      }
    p[i,num.sources] <- 1-sum(p[i,1:(num.sources-1)]);
   
    ##rescale p.pop for concentration dependence
    p.inddenom[i,1:num.iso] <- mu.conc2%*%p[i,]
    for(iso in 1:num.iso) {
      for(source in 1:num.sources) {
        pIso.ind[i,iso,source]  <- mu.conc2[iso,source]*p[i,source]/p.inddenom[i,iso]
      }
    }
 }#end for i
  
  #####################
  ##Source Estimation##
  #####################
  ##estimate sources
  for(source in 1:num.sources) { 
   ##covariance matrix
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
			tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
		tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(source2 in 1:num.iso) {
		tau.source.temp[source2,source2,source] ~ dgamma(.0010,.0010)
    }
       
    ##build source correlation matrix
   for(sourcex in 1:num.iso) {
      for(sourcey in 1:num.iso) {
 		  rho.source[sourcex,sourcey,source] ~ dunif(-0.99,0.99);
		}
    }
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        rho.mat[source, sourcex,sourcey] <- rho.source[sourcex, sourcey, source]*rho.flag
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        rho.mat[source,sourcex,sourcey] <- rho.mat[source, sourcey, sourcex];
      }
    }
    for(source2 in 1:num.iso) {
      rho.mat[source,source2,source2] <- 1;
    }      

    cov.source[1:num.iso,1:num.iso,source] <- inverse(tau.source.temp[,,source])
    #tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[source,,]%*%cov.source[,,source]))
    tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[source,,]%*%cov.source[,,source]) )
 
    ##draws subsource means and fit to data
    for(sub in 1:subsource.vec[source]) { 
		subsource[source, sub, 1:num.iso] ~ dmnorm(mu.prior.mu, mu.prior.cov)
	}
	
	for(iso in 1:num.iso) {
		mu.source[source, iso] <- mean( subsource[source, 1:subsource.vec[source], iso])
	}
		
	for(sub in 1:subsource.vec[source]) {
		for(index  in subsource.samples[source,sub,1]:subsource.samples[source,sub,2]) {
			source.mat[index, ] ~ dmnorm(subsource[source,sub, ], tau.source[ , ,source]);
		}  
     }
  }#end for sources
  
  
  ####################### 
  ##draw residual error##
  #######################
  for(i in 1:(num.iso-1)) {
	  for(j in (i+1):num.iso) {
		res.tau[i,j] <- 0;
		res.tau[j,i] <- 0;
		}
	}
  for(iso in 1:num.iso) {
    res.tau[iso,iso] ~ dgamma(1e-3,1e-3)
  }
  res.err[1:num.iso,1:num.iso] <- inverse(res.tau)

  ##rescale sources by p
  for(i in 1:N) {
    ##rescale covariance and include fractionation
    for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
      #  covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source])*pIso.ind[i,iso,source]
        covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source] )*pIso.ind[i,iso,source]
      }
      for(isox in 2:(num.iso)) {
    for(isoy in 1:(isox-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
      for(isoy in 2:(num.iso)) {
    for(isox in 1:(isoy-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
  
      obscov.mat[1:num.iso,1:num.iso,source,i] <- (covfrac.source[,,source,i]%*%rho.mat[source,,]%*%covfrac.source[,,source,i] +  res.err)  
  
    }#end for source

    for(x in 1:num.iso) {
      for(y in 1:num.iso) {
		sumobscov.mat[x,y,i] <- sum(obscov.mat[x,y,1:num.sources,i])  
      }
    }
    for(iso in 1:num.iso) {
      mu.mix[i,iso] <- pIso.ind[i,iso,]%*%(mu.source[,iso])
    }
  }#end for i

  ##get the sd's for jack
  for(iso in 1:num.iso) {
  
    sd.res[iso] <- sqrt(res.err[iso,iso])
            
    for(source in 1:num.sources) {

        sd.source[source,iso] <- sqrt(cov.source[iso,iso,source])
        sd.conc[source,iso] <- 1/sqrt(D.tau[iso,iso,source])

    }
  }

  mu.conc	<- t(mu.conc2)

  ##calculate the likelihoods for the N individuals.
  for(ind in 1:N) {
    mix.prcsn[1:num.iso,1:num.iso,ind] <- inverse(sumobscov.mat[,,ind]  )
    for(j in 1:ind.counts[ind]) {      
      ind.array[1:num.iso,ind,j] ~ dmnorm(mu.mix[ind,1:num.iso], mix.prcsn[1:num.iso,1:num.iso, ind]);    
    }
  }

}#end model
"

