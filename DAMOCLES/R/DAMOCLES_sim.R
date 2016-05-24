roundn = function(x, digits = 0)
{
    fac = 10^digits
    return(trunc(fac * x + 0.5)/fac)
}

sample2 = function(x,size,replace = FALSE,prob = NULL)
{
    if(length(x) == 1)
    { 
        x = c(x,x)
        prob = c(prob,prob)
    }
    return(sample(x,size,replace,prob))
}

DAMOCLES_sim = function(
	phy,
	gamma_0,
	gamma_td,
	mu,
	sigma,
	psiBranch,
	psiTrait,
	z,
	phi,
	traitOpt,
	br0,
	br_td,
	nTdim,
	root.state,
	root.trait.state,
	plotit = FALSE,
	keepExtinct = FALSE
	)
{
	  
  dn = dist.nodes(phy) # find distances between nodes
	ntips = length(phy$tip.label) # number of tips
	nbranch = 2 * ntips - 2 # number of branches
	# create a dataframe to store presence/absence values of each edge in the tree
  patable = data.frame(p = phy$edge[,1],d = phy$edge[,2],age.start = rep(NA,nbranch),age.end = rep(NA,nbranch),extant = rep(0,nbranch),state = rep(NA,nbranch))
	patable$age.start = dn[patable$p,ntips + 1] # times when each edge starts
	patable$age.end = dn[patable$d,ntips + 1] # times when each edge ends

	patable$age.end[which(patable$d < length(phy$tip.label) + 1)] = max(patable$age.end)
		
	ce = which(patable$p == min(patable$p))# find the edges diverging from the crown
	patable$extant[ce] = 1  # initialise these edges as extant
	patable$state = rep(NA,dim(patable)[1]) # initialise the states i.e. presence/absence of extant edges
	patable$state[ce] = sample(c(0,root.state))
	
	patable$tips[which(patable$d <= ntips)] = 1
  patable$y[which(patable$tips == 1)] = patable$d[which(patable$tips == 1)]
    
  if(plotit == TRUE)
  {
     plot.phylo(phy,main = "present (red), absent (blue)")
		 while(length(na.omit(patable$y)) < length(patable$y))
     {
	    	for(i in 1:length(patable$y))
        {
		    	 if(!is.na(patable$y[i]))
           {
			     	  focalp = patable$p[i]
			       	focaly = patable$y[which(patable$p == focalp)]
  		     	  if(length(na.omit(focaly)) == 2)
              {
	    	    	   focald = which(patable$d == focalp)
			    	 	   patable$y[focald] = (focaly[1] + focaly[2])/2
   	   		    }
			     }
			  }
     }
  }	
	
	if(nTdim > 1)
  { # store traits as a vector or a matrix for 1 or multiple dimensions
		traits = matrix(NA,ncol = dim(patable)[1],nrow = nTdim)
		traits[,ce] = 0
	} else {
		traits = rep(NA,dim(patable)[1])
		traits[ce] = 0	
	}
		
	tstep = min(patable$age.start)
  tlastevent = max(patable$age.end)
		
	while(tstep < tlastevent)
  {		
		extant = which(patable$extant == 1) # find extant edges
		extant.id = patable$d[extant]
		pa = patable$state[extant] # extract the states of extant edges
		
		Itot.max = gamma_0 * length(pa) # use maximum rates to calculate waiting time till next event
		Etot.max = mu * length(pa)
		
		tstep0 = tstep
		wt = -log(runif(1))/(Itot.max + Etot.max)		# waiting time untill the next immigration/extinction event
		tstep1 = tstep0 + wt		# calculate the time when the immigration/extinction event will occur
    tnextevent = min(patable$age.end[extant])
    tstep = min(tstep1,tnextevent)
		#tstep = min(tstep1, min(patable$age.end[extant]))		# if it occurs after the next branching or extinction event then set tstep to this time
			
		if(br_td == 0)
    {
			wt = tstep - tstep0 # calculate the realised waiting time
		} else {
			wt = (max(patable$age.end) * (1 - exp(br_td*tstep))) - (max(patable$age.end) * (1 - exp(br_td * tstep0)))
		}	
		
		if(plotit == TRUE)
    {
	     for(i in unique(extant))
       {
	     		if(patable$state[i] == 1)
           {
		       		lines(c(tstep0,tstep),c(patable$y[i],patable$y[i]),col = "red",lwd = 2)
           } else {
	            lines(c(tstep0,tstep),c(patable$y[i],patable$y[i]),col = "blue",lwd = 2)
       	   }	
       }     
    }
				
		sdnorm = sqrt(br0 * wt)
		traits = traits + rnorm(length(traits),mean = 0,sd = sdnorm) # simulate brownian trait evolution
				
		if(tstep < tlastevent)
    { # if we have yet to reach the present then simulate the next event 
			if(tstep < tnextevent)
      { # if the next event is an immigration or local extinction event 
				
				Itot.real = 0
				
				if(length(pa) - sum(pa) > 0)
        { # if there are species absent from the community
					
					Dpaynotpa = 1 - pa
					Itot.real = gamma_0 * sum(Dpaynotpa)
					
					if(gamma_td > 0)
          { # if we have a slowdown in immigration rates.......
						Dpaynotpa = 1 - pa
						Dpaynotpa = gamma_0/(1 + (gamma_td * tstep)) * Dpaynotpa
						Itot.real = sum(Dpaynotpa)
					}
					
					if(psiBranch > 0 & sum(pa) > 0)
          { #	if we have phylogenetic repulsion and species are present in the community already...........
						Db = matrix(rep(patable$age.end[extant] - tstep, sum(patable$extant)), nrow = sum(patable$extant))#find the distance between the current time step and the end of each extant edge - create a matrix
						futD = Db + t(Db)#transpose and add to the original matrix giving us the combined distance from the current time step to the ends of each pair of edges #i.e. call it futD (future distance)
						D = dn[extant.id, extant.id] - futD #the distance between two branches at the current time step is then the distance between the end of the edges minus the future distance
						D = replace(D, D < 0, 0)  # due to rounding error we can get negative distances so set these to zero 
						D[,which(pa == 1)] = NaN
            D[which(pa == 0),] = NaN
						Dpay = D^z/(psiBranch^z + D^z)
						Dpaynotpa = colProds(Dpay,na.rm = TRUE) * (1 - pa)
						Dpaynotpa = gamma_0 * Dpaynotpa
						Itot.real = sum(Dpaynotpa)
					} 
										
					if(phi > 0)
          { # if we have environmental filtering.......
						if(nTdim > 1)
            { # calculate distance of species from trait optima
							D = traitOpt - traits[,extant]
							D = sqrt(colSums(D^2))
						} else {	
							D = abs(traitOpt - traits[extant])
						}
						Dpaynotpa = dexp(D,rate = phi)/dexp(0,rate = phi)
						Dpaynotpa = gamma_0 * Dpaynotpa
						Dpaynotpa = Dpaynotpa * (1 - pa)
						Itot.real = sum(Dpaynotpa)
					}	
					
					if(psiTrait > 0 & sum(pa) > 0)
          {  # if we have trait repulsion.......					
						if(nTdim > 1)
            {
							D = dist(t(traits[,extant]), method = "euclidean", diag = FALSE, upper = TRUE, p = 2)
						}else{
							D = dist(traits[extant], method = "euclidean", diag = FALSE, upper = TRUE, p = 2)
						}	
						D = as.matrix(D)
							
						D[which(pa == 0),] = NaN
						Dpay = D^z/(psiTrait^z + D^z)
						Dpaynotpa  = colProds(Dpay,na.rm = TRUE) * (1 - pa)
						Dpaynotpa = gamma_0 * Dpaynotpa
						Itot.real = sum(Dpaynotpa)							
					}
				}
				
				Etot.real = mu * sum(pa)
				Itot = Itot.real /(Itot.max + Etot.max) # probabilities of events are rates at the time step divided by the total maximum rate
				Etot = Etot.real /(Itot.max + Etot.max)
					
				event = sample(x = 0:2, size = 1, replace = F, prob = c(1 - (Itot + Etot), Itot, Etot)) # then choose which event occurs 

				# simulate the next event - immigration, local extinction, speciation, global extinction or nothing

				if(event == 1)
        { # if it's an immigration event.....
					patable$state[sample2(extant, 1, replace = FALSE, prob = Dpaynotpa)] = 1					
				} 
				if(event == 2)
        {	# if it's a local extinction event........
          patable$state[sample2(extant[which(pa  == 1)], 1)] = 0					
				} 			
			} else {	# if it's a speciation or global extinction event....... 
				focedge = extant[which(patable$age.end[extant] == min(patable$age.end[extant]))] # find which edge this event occurs along	
				for(fe in unique(focedge))
        { # there may be multiple events at the same time  - so loop through these
					dedge = which(patable$p == patable$d[fe])
					if(length(dedge) > 0)
          {	# if it's a speciation event......
						patable$extant[fe] = 0		# then assign this edge as extinct
						patable$extant[dedge] = 1		# and assign daughter edges as extant
						if(nTdim > 1)
            {
							traits[,dedge] = traits[,fe] # give the daughters the trait of the parent
						} else {
							traits[dedge] = traits[fe]
						}	
						ds = sample(c(0,patable$state[fe]),1,prob = c(1 - sigma,sigma)) # assign the speciation event as in-situ or not given sigma
						patable$state[dedge] = sample(c(ds,patable$state[fe]))	# then assign geographic states at random (needed if the speciation event is allopatric) 
					} else {	# if it is an extinction event......
						patable$extant[fe] = 0	# then assign this edge as extinct
					}
				}
			}
		}
	}	
	patablelist = list()
	if(keepExtinct == FALSE)
  {
		if(nTdim > 1)
    { # prune out extinct species
			traits  = traits[,which(patable$extant == 1)]
		} else {
			traits  = traits[which(patable$extant == 1)]
		}
		 
		patable  = patable[which(patable$extant == 1),] # prune out extinct species
	}
	patablelist[[1]] = patable
	patablelist[[2]] = traits
	return(patablelist)	
}					

					