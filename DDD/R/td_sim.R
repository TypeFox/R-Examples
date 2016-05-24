td_sim = function(pars,age,ddmodel = 1,methode = 'ode45')
  # Other methdos: 'ode45', 'lsodes'
{
  # Simulation of diversity-dependent process
  #  . start from crown age
  #  . no additional species at crown node
  #  . no missing species in present
  # pars = [la mu K]
  #  . la = speciation rate per species
  #  . mu = extinction rate per species
  #  . K = diversification carrying capacity
  #  . r = ratio of diversity-dependence in extinction rate over speciation rate
  # age = crown age
  # ddmodel = mode of diversity-dependence
  #  . ddmodel == 1 : linear dependence in speciation rate with parameter K
  #  . ddmodel == 1.3: linear dependence in speciation rate with parameter K'
  #  . ddmodel == 2 : exponential dependence in speciation rate
  #  . ddmodel == 2.1: variant with offset at infinity
  #  . ddmodel == 2.2: 1/n dependence in speciation rate
  #  . ddmodel == 2.3: exponential dependence in speciation rate with parameter x
  #  . ddmodel == 3 : linear dependence in extinction rate
  #  . ddmodel == 4 : exponential dependence in extinction rate
  #  . ddmodel == 4.1: variant with offset at infinity
  #  . ddmodel == 4.2: 1/n dependence in speciation rate
  #  . ddmodel == 5 : linear dependence in speciation rate and in extinction rate
  
  # initialisation of the time dependent function
  
  # pars1 contains model parameters
  # - pars1[1] = la0 = speciation rate
  # - pars1[2] = mu0 = extinction rate
  # - pars1[3] = la1 = parameter in exponential decay of speciation rate, or K in diversity-dependence-like models (default 0)
  # - pars1[4] = mu2 = parameter in exponential decay of extinction rate (default 0)
  # - pars1[5] = T0 = age at which lambda is lambda0 (default T0 = age of phylogeny)
  
  la0 = pars[1]
  mu0 = pars[2]
  K = pars[3]
  
  #test on the parameters
  
  if(la0 < mu0) {stop('the function is designed for lambda_0 > mu_0')}
  if(mu0 < 0) {stop('per species rates should be positive')}
  if(K < 1) {stop('clade level carrying capacity should be positive')}
  
  Kprime = la0/(la0 - mu0) * K
  
  soc = 2
  nbeq = 1000
  tdmodel = 4
   
  done = 0
  while(done == 0)
  {
    # number of species N at time t
    # i = index running through t and N
    t = rep(0,1)  
    L = matrix(0,2,4)
        
    i = 1 # counter variable for the potential events    
    ev = 1 #counter variable for the events events (conditioning is at crown age)
    t[1] = 0    #initialisation of a time vector for all potential events
    N = 2
    
    lx = min(1 + ceiling(Kprime),nbeq)
    variables = rep(0,lx)
    variables[soc + 1] = 1
    
    # lambda(t) at t=0, time dependent model
    latd = mu0 + ((la0 - mu0) * soc - soc^2*(la0-mu0)/K)/soc 
    mutd = mu0
    
    #firt upper values for the rates
    lamax = latd 
    mumax = mu0             
    # L = data structure for lineages,
    # . L[,1] = branching times
    # . L[,2] = index of parent species
    # . L[,3] = index of daughter species
    # . L[,4] = time of extinction
    # j = index running through L
    L[1,1:4] = c(0,0,-1,-1)
    L[2,1:4] = c(0,-1,2,-1)
    linlist = c(-1,2)
    newL = 2
    
    # parameter next potential event
    denom = (lamax + mumax) * N[i]
    
    #update the time with the first potential event
    t[i + 1] = t[i] + rexp(1,denom)
    
    while(t[i + 1] <= age)
    {
      i = i + 1     
      
      # sample an individual 
      ranL = sample(linlist,1)
      
      # calculation of the rate latd at the new time step
      y = ode(variables, c(t[i-1], t[i]), td_loglik_rhs_sim,c(pars[1:min(4,length(pars))],tdmodel-3, lx), rtol = 1e-10, atol = 1e-16, method = methode)
      
      variables = y[2, 2:(lx + 1)]
      expn = sum((0:(lx - 1)) * variables[1:lx])
      expn2 = sum((0:(lx - 1))^2 * variables[1:lx])
      dEN_dt = (la0 - mu0) * expn - expn2 * (la0 - mu0)/K
      
      latd = mu0 + dEN_dt/expn
            
      if(latd > lamax)
      {
         stop('latd should be a decreasing function of time')
      }      
      if( ((latd + mutd)/(lamax + mumax)) >= runif(1) )  #does the next potential event occur?                    
      {
        ev = ev + 1
        lamax = latd
        
        if( latd/(latd + mutd) >= runif(1) )  # Is it a speciation event?
        {           
          N[ev] = N[ev - 1] + 1
          newL = newL + 1
          L = rbind(L,c(t[i],ranL,sign(ranL) * newL,-1))  
          linlist = c(linlist,sign(ranL) * newL)
        } else {         
          N[ev] = N[ev - 1] - 1
          L[abs(ranL),4] = t[i]
          w = which(linlist == ranL)
          linlist = linlist[-w]
          linlist = sort(linlist)          
        }
      } 
   
      # Is one of the two crown lineages extinct ?
      if(sum(linlist < 0) == 0 | sum(linlist > 0) == 0)
      {
        t[i + 1] = Inf
      } 
      else 
      {
        #parameter for the next potential event is updated with new lambda_max, mu_max, and number of species      
        denom = (lamax + mumax) * N[ev]
        #time is updated by adding a new potential step
        t[i + 1] = t[i] + rexp(1,denom)
      } 
    }
    
    if(sum(linlist < 0) == 0 | sum(linlist > 0) == 0)
    {
      done = 0
    } else {
      done = 1
    }
  }
  
  L[,1] = age - c(L[,1])
  notmin1 = which(L[,4] != -1)
  L[notmin1,4] = age - c(L[notmin1,4])
  L[which(L[,4] == age + 1),4] = -1
  tes = L2phylo(L,dropextinct = T)
  tas = L2phylo(L,dropextinct = F)
  out = list(tes,tas,L)
  return(out)  
}